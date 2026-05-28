//! Rack, arrangement, gel-preview, and rack-label GUI helpers.
//!
//! This module is a move-only extraction from `app.rs`: it keeps the rack
//! workspace, arrangement gel preview, rack-label preview, and placement
//! dialogs close to `GENtleApp` while reducing the top-level app monolith.

use super::*;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum RackDragState {
    Sample {
        from_coordinate: String,
        arrangement_id: String,
        role_label: String,
    },
    Samples {
        arrangement_id: String,
        from_coordinates: Vec<String>,
    },
    ArrangementBlock {
        arrangement_id: String,
        from_coordinate: String,
    },
    ArrangementBlocks {
        arrangement_ids: Vec<String>,
        from_coordinate: String,
    },
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct RackGhostPreviewCell {
    pub(super) predicted_entry: Option<crate::engine::RackPlacementEntry>,
}

#[derive(Clone)]
pub(super) struct RackDropGhostState {
    pub(super) rack_id: String,
    pub(super) changed_coordinates: BTreeSet<String>,
    pub(super) started_at: Instant,
}

#[derive(Clone)]
pub(super) struct ArrangementGelPreviewState {
    pub(super) arrangement_id: String,
    pub(super) arrangement_title: String,
    pub(super) saved_ladders: Vec<String>,
    pub(super) left_ladder_name: String,
    pub(super) right_ladder_name: String,
    pub(super) agarose_percent_text: String,
    pub(super) buffer_model: crate::engine::GelBufferModel,
    pub(super) topology_aware: bool,
    pub(super) svg_uri: String,
    pub(super) status: String,
}

impl Default for ArrangementGelPreviewState {
    fn default() -> Self {
        let defaults = crate::engine::GelRunConditions::default();
        Self {
            arrangement_id: String::new(),
            arrangement_title: String::new(),
            saved_ladders: vec![],
            left_ladder_name: String::new(),
            right_ladder_name: String::new(),
            agarose_percent_text: format!("{:.1}", defaults.agarose_percent),
            buffer_model: defaults.buffer_model,
            topology_aware: defaults.topology_aware,
            svg_uri: String::new(),
            status: String::new(),
        }
    }
}

#[derive(Clone, Default)]
pub(super) struct RackLabelsPreviewState {
    pub(super) rack_id: String,
    pub(super) rack_title: String,
    pub(super) arrangement_id: Option<String>,
    pub(super) arrangement_title: Option<String>,
    pub(super) svg_uri: String,
    pub(super) status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub(super) struct PersistedRackWorkspace {
    pub(super) help_strip_collapsed: bool,
    pub(super) help_strip_pinned_open: bool,
    pub(super) help_strip_successful_move_count: u32,
    pub(super) help_strip_auto_minimized: bool,
}

impl GENtleApp {
    pub(super) fn prompt_export_serial_gel_svg(
        &mut self,
        default_stem: &str,
        container_ids: Option<Vec<String>>,
        arrangement_id: Option<String>,
        ladders: Option<Vec<String>>,
        conditions: Option<crate::engine::GelRunConditions>,
    ) {
        let stem = Self::sanitize_file_stem(default_stem, "serial_gel");
        let default_file_name = format!("{stem}.gel.svg");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Serial gel SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::RenderPoolGelSvg {
                inputs: vec![],
                path: path_text.clone(),
                ladders,
                container_ids,
                arrangement_id,
                conditions,
            });
        match result {
            Ok(op_result) => {
                self.app_status = op_result
                    .messages
                    .first()
                    .cloned()
                    .unwrap_or_else(|| format!("Wrote serial gel SVG to '{path_text}'"));
            }
            Err(e) => {
                self.app_status = format!("Could not export serial gel SVG: {}", e.message);
            }
        }
    }

    pub(super) fn arrangement_gel_preview_title(&self) -> String {
        if self
            .arrangement_gel_preview
            .arrangement_title
            .trim()
            .is_empty()
        {
            "Arrangement Gel".to_string()
        } else {
            format!(
                "Arrangement Gel — {}",
                self.arrangement_gel_preview.arrangement_title.trim()
            )
        }
    }

    pub(super) fn describe_arrangement_ladders(ladders: &[String]) -> String {
        if ladders.is_empty() {
            "auto".to_string()
        } else {
            ladders.join(" + ")
        }
    }

    pub(super) fn set_arrangement_gel_preview_controls_from_ladders(&mut self, ladders: &[String]) {
        self.arrangement_gel_preview.saved_ladders = ladders.to_vec();
        self.arrangement_gel_preview.left_ladder_name.clear();
        self.arrangement_gel_preview.right_ladder_name.clear();
        match ladders {
            [] => {}
            [name] => {
                self.arrangement_gel_preview.left_ladder_name = name.clone();
                self.arrangement_gel_preview.right_ladder_name = name.clone();
            }
            [left, right, ..] => {
                self.arrangement_gel_preview.left_ladder_name = left.clone();
                self.arrangement_gel_preview.right_ladder_name = right.clone();
            }
        }
    }

    pub(super) fn arrangement_gel_preview_effective_ladders(&self) -> Vec<String> {
        let left = self.arrangement_gel_preview.left_ladder_name.trim();
        let right = self.arrangement_gel_preview.right_ladder_name.trim();
        match (left.is_empty(), right.is_empty()) {
            (true, true) => vec![],
            (false, true) => vec![left.to_string()],
            (true, false) => vec![right.to_string()],
            (false, false) if left == right => vec![left.to_string()],
            (false, false) => vec![left.to_string(), right.to_string()],
        }
    }

    pub(super) fn arrangement_gel_preview_effective_conditions(
        &self,
    ) -> std::result::Result<crate::engine::GelRunConditions, String> {
        let raw = self.arrangement_gel_preview.agarose_percent_text.trim();
        let agarose_percent = if raw.is_empty() {
            crate::engine::GelRunConditions::default().agarose_percent
        } else {
            raw.parse::<f32>()
                .map_err(|e| format!("Invalid agarose percent '{raw}': {e}"))?
        };
        Ok(crate::engine::GelRunConditions {
            agarose_percent,
            buffer_model: self.arrangement_gel_preview.buffer_model,
            topology_aware: self.arrangement_gel_preview.topology_aware,
        }
        .normalized())
    }

    pub(super) fn coerce_arrangement_gel_preview_pair(&mut self) {
        let left = self
            .arrangement_gel_preview
            .left_ladder_name
            .trim()
            .to_string();
        let right = self
            .arrangement_gel_preview
            .right_ladder_name
            .trim()
            .to_string();
        if left.is_empty() && !right.is_empty() {
            self.arrangement_gel_preview.left_ladder_name = right.clone();
        } else if !left.is_empty() && right.is_empty() {
            self.arrangement_gel_preview.right_ladder_name = left.clone();
        }
    }

    pub(super) fn open_arrangement_gel_preview_dialog(&mut self, arrangement_id: &str) {
        let arrangement_id = arrangement_id.trim();
        if arrangement_id.is_empty() {
            self.app_status =
                "Arrangement gel preview requires a non-empty arrangement id".to_string();
            return;
        }
        if self.show_arrangement_gel_preview_dialog
            && self.arrangement_gel_preview.arrangement_id == arrangement_id
        {
            self.mark_window_open_or_focus(Self::arrangement_gel_preview_viewport_id(), true);
            return;
        }
        self.refresh_lineage_cache_if_needed();
        let Some(arrangement) = self
            .lineage_arrangements
            .iter()
            .find(|row| row.arrangement_id == arrangement_id)
            .cloned()
        else {
            self.app_status = format!("Arrangement '{arrangement_id}' is not available");
            return;
        };
        self.arrangement_gel_preview.arrangement_id = arrangement.arrangement_id.clone();
        self.arrangement_gel_preview.arrangement_title = if arrangement.name.trim().is_empty() {
            arrangement.arrangement_id.clone()
        } else {
            arrangement.name.clone()
        };
        let default_conditions = crate::engine::GelRunConditions::default();
        self.arrangement_gel_preview.agarose_percent_text =
            format!("{:.1}", default_conditions.agarose_percent);
        self.arrangement_gel_preview.buffer_model = default_conditions.buffer_model;
        self.arrangement_gel_preview.topology_aware = default_conditions.topology_aware;
        self.arrangement_gel_preview.status.clear();
        self.arrangement_gel_preview.svg_uri.clear();
        self.set_arrangement_gel_preview_controls_from_ladders(&arrangement.ladders);
        if arrangement.ladders.len() > 2 {
            self.arrangement_gel_preview.status = format!(
                "Saved arrangement stores {} ladder names; this preview editor currently exposes left/right flanks only.",
                arrangement.ladders.len()
            );
        }
        self.show_arrangement_gel_preview_dialog = true;
        self.mark_window_open_or_focus(Self::arrangement_gel_preview_viewport_id(), false);
        self.refresh_arrangement_gel_preview_svg();
    }

    pub(super) fn refresh_arrangement_gel_preview_svg(&mut self) {
        let arrangement_id = self
            .arrangement_gel_preview
            .arrangement_id
            .trim()
            .to_string();
        if arrangement_id.is_empty() {
            self.arrangement_gel_preview.status =
                "No arrangement is selected for gel preview.".to_string();
            self.arrangement_gel_preview.svg_uri.clear();
            return;
        }
        let override_ladders = self.arrangement_gel_preview_effective_ladders();
        let conditions = match self.arrangement_gel_preview_effective_conditions() {
            Ok(conditions) => conditions,
            Err(err) => {
                self.arrangement_gel_preview.status = err;
                self.arrangement_gel_preview.svg_uri.clear();
                return;
            }
        };
        let layout_result = {
            let engine = self.engine.read().unwrap();
            engine.build_serial_gel_layout_for_render(
                &[],
                None,
                Some(arrangement_id.as_str()),
                Some(override_ladders.as_slice()),
                Some(&conditions),
            )
        };
        let layout = match layout_result {
            Ok(layout) => layout,
            Err(err) => {
                self.arrangement_gel_preview.status =
                    format!("Could not refresh arrangement gel preview: {}", err.message);
                self.arrangement_gel_preview.svg_uri.clear();
                return;
            }
        };
        let svg = crate::pool_gel::export_pool_gel_svg(&layout);
        let stamp = now_unix_ms_u64();
        let png_path = env::temp_dir().join(format!("gentle_arrangement_gel_preview_{stamp}.png"));
        let raster_result = (|| -> std::result::Result<(), String> {
            let mut opt = usvg::Options {
                resources_dir: png_path.parent().map(|path| path.to_path_buf()),
                ..usvg::Options::default()
            };
            opt.fontdb_mut().load_system_fonts();
            let tree = usvg::Tree::from_str(&svg, &opt)
                .map_err(|e| format!("Could not parse arrangement gel SVG: {e}"))?;
            let pixmap_size = tree
                .size()
                .to_int_size()
                .scale_by(1.0)
                .ok_or_else(|| "Could not determine arrangement gel raster size".to_string())?;
            let mut pixmap = tiny_skia::Pixmap::new(pixmap_size.width(), pixmap_size.height())
                .ok_or_else(|| {
                    format!(
                        "Could not allocate arrangement gel raster {}x{}",
                        pixmap_size.width(),
                        pixmap_size.height()
                    )
                })?;
            resvg::render(&tree, tiny_skia::Transform::default(), &mut pixmap.as_mut());
            pixmap.save_png(&png_path).map_err(|e| {
                format!(
                    "Could not write temporary arrangement gel preview PNG '{}': {e}",
                    png_path.display()
                )
            })?;
            Ok(())
        })();
        match raster_result {
            Ok(()) => {
                self.arrangement_gel_preview.svg_uri = Self::file_uri_from_path(&png_path);
                self.arrangement_gel_preview.status = format!(
                    "Previewing {} sample lane(s), {} sequence(s); ladders: {}",
                    layout.sample_count,
                    layout.pool_member_count,
                    Self::describe_arrangement_ladders(&layout.selected_ladders)
                );
            }
            Err(err) => {
                self.arrangement_gel_preview.status = err;
                self.arrangement_gel_preview.svg_uri.clear();
            }
        }
    }

    pub(super) fn save_arrangement_gel_preview_ladders(&mut self) {
        let arrangement_id = self
            .arrangement_gel_preview
            .arrangement_id
            .trim()
            .to_string();
        if arrangement_id.is_empty() {
            self.arrangement_gel_preview.status =
                "No arrangement is selected for ladder updates.".to_string();
            return;
        }
        let ladders = self.arrangement_gel_preview_effective_ladders();
        let op_result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::SetArrangementLadders {
                arrangement_id: arrangement_id.clone(),
                ladders: (!ladders.is_empty()).then_some(ladders.clone()),
            });
        match op_result {
            Ok(result) => {
                self.app_status =
                    result.messages.first().cloned().unwrap_or_else(|| {
                        format!("Updated arrangement '{arrangement_id}' ladders")
                    });
                self.arrangement_gel_preview.status = self.app_status.clone();
                self.set_arrangement_gel_preview_controls_from_ladders(&ladders);
                self.lineage_cache_valid = false;
                self.refresh_lineage_cache_if_needed();
                self.refresh_arrangement_gel_preview_svg();
            }
            Err(err) => {
                self.arrangement_gel_preview.status =
                    format!("Could not update arrangement ladders: {}", err.message);
            }
        }
    }

    pub(super) fn rack_labels_preview_title(&self) -> String {
        let rack_label = if self.rack_labels_preview.rack_title.trim().is_empty() {
            self.rack_labels_preview.rack_id.trim().to_string()
        } else {
            self.rack_labels_preview.rack_title.trim().to_string()
        };
        match self
            .rack_labels_preview
            .arrangement_title
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            Some(arrangement_title) => {
                format!("Rack Labels — {} / {}", rack_label, arrangement_title)
            }
            None if rack_label.is_empty() => "Rack Labels".to_string(),
            None => format!("Rack Labels — {}", rack_label),
        }
    }

    pub(super) fn open_rack_labels_preview_dialog(
        &mut self,
        rack_id: &str,
        arrangement_id: Option<&str>,
    ) {
        let rack_id = rack_id.trim();
        if rack_id.is_empty() {
            self.app_status = "Rack label preview requires a non-empty rack id".to_string();
            return;
        }
        self.refresh_lineage_cache_if_needed();
        let Some(rack) = self
            .lineage_racks
            .iter()
            .find(|row| row.rack_id == rack_id)
            .cloned()
        else {
            self.app_status = format!("Rack '{}' is not available", rack_id);
            return;
        };
        let arrangement_id = arrangement_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let arrangement_title = arrangement_id.as_ref().map(|id| {
            self.lineage_arrangements
                .iter()
                .find(|row| row.arrangement_id == *id)
                .map(|row| {
                    if row.name.trim().is_empty() {
                        row.arrangement_id.clone()
                    } else {
                        row.name.clone()
                    }
                })
                .unwrap_or_else(|| id.clone())
        });
        if self.show_rack_labels_preview_dialog
            && self.rack_labels_preview.rack_id == rack_id
            && self.rack_labels_preview.arrangement_id == arrangement_id
        {
            self.mark_window_open_or_focus(Self::rack_labels_preview_viewport_id(), true);
            return;
        }
        self.rack_labels_preview.rack_id = rack.rack_id.clone();
        self.rack_labels_preview.rack_title = if rack.name.trim().is_empty() {
            rack.rack_id.clone()
        } else {
            format!("{} ({})", rack.name.trim(), rack.rack_id.trim())
        };
        self.rack_labels_preview.arrangement_id = arrangement_id;
        self.rack_labels_preview.arrangement_title = arrangement_title;
        self.rack_labels_preview.status.clear();
        self.rack_labels_preview.svg_uri.clear();
        self.show_rack_labels_preview_dialog = true;
        self.mark_window_open_or_focus(Self::rack_labels_preview_viewport_id(), false);
        self.refresh_rack_labels_preview_svg();
    }

    pub(super) fn refresh_rack_labels_preview_svg(&mut self) {
        let rack_id = self.rack_labels_preview.rack_id.trim().to_string();
        if rack_id.is_empty() {
            self.rack_labels_preview.status = "No rack is selected for label preview.".to_string();
            self.rack_labels_preview.svg_uri.clear();
            return;
        }
        let arrangement_id = self
            .rack_labels_preview
            .arrangement_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let rendered = {
            let engine = self.engine.read().unwrap();
            engine.render_rack_labels_svg_text(
                &rack_id,
                arrangement_id.as_deref(),
                self.rack_label_sheet_preset,
            )
        };
        let (svg, label_count) = match rendered {
            Ok(result) => result,
            Err(err) => {
                self.rack_labels_preview.status =
                    format!("Could not refresh rack label preview: {}", err.message);
                self.rack_labels_preview.svg_uri.clear();
                return;
            }
        };
        let scope_stem = arrangement_id
            .as_deref()
            .map(|value| Self::sanitize_file_stem(value, "arrangement"))
            .unwrap_or_else(|| "whole_rack".to_string());
        let base_name = format!(
            "gentle_rack_labels_preview_{}_{}_{}",
            Self::sanitize_file_stem(&rack_id, "rack"),
            scope_stem,
            self.rack_label_sheet_preset.as_str()
        );
        let svg_path = std::env::temp_dir().join(format!("{base_name}.svg"));
        let png_path = std::env::temp_dir().join(format!("{base_name}.png"));
        let raster_result = (|| -> std::result::Result<(), String> {
            fs::write(&svg_path, &svg).map_err(|e| {
                format!(
                    "Could not write temporary rack labels preview SVG '{}': {e}",
                    svg_path.display()
                )
            })?;
            let mut opt = usvg::Options::default();
            opt.fontdb_mut().load_system_fonts();
            let tree = usvg::Tree::from_str(&svg, &opt)
                .map_err(|e| format!("Could not parse rack labels SVG: {e}"))?;
            let pixmap_size = tree
                .size()
                .to_int_size()
                .scale_by(1.0)
                .ok_or_else(|| "Could not determine rack label raster size".to_string())?;
            let mut pixmap = tiny_skia::Pixmap::new(pixmap_size.width(), pixmap_size.height())
                .ok_or_else(|| {
                    format!(
                        "Could not allocate rack label raster {}x{}",
                        pixmap_size.width(),
                        pixmap_size.height()
                    )
                })?;
            resvg::render(&tree, tiny_skia::Transform::default(), &mut pixmap.as_mut());
            pixmap.save_png(&png_path).map_err(|e| {
                format!(
                    "Could not write temporary rack labels preview PNG '{}': {e}",
                    png_path.display()
                )
            })?;
            Ok(())
        })();
        match raster_result {
            Ok(()) => {
                self.rack_labels_preview.svg_uri = Self::file_uri_from_path(&png_path);
                let scope_label = self
                    .rack_labels_preview
                    .arrangement_title
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(|value| format!(" arrangement '{}'", value))
                    .unwrap_or_else(|| " whole rack scope".to_string());
                self.rack_labels_preview.status = format!(
                    "Previewing {} rack label card(s) for{} using preset '{}'.",
                    label_count,
                    scope_label,
                    self.rack_label_sheet_preset.as_str()
                );
            }
            Err(err) => {
                self.rack_labels_preview.status = err;
                self.rack_labels_preview.svg_uri.clear();
            }
        }
    }

    pub(super) fn rack_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("gentle_rack_view")
    }

    pub(super) fn rack_profile_label(kind: RackProfileKind) -> &'static str {
        match kind {
            RackProfileKind::SmallTube4x6 => "Small tube rack (4 x 6)",
            RackProfileKind::Plate6 => "Cell culture plate 6-well",
            RackProfileKind::Plate12 => "Cell culture plate 12-well",
            RackProfileKind::Plate24 => "Cell culture plate 24-well",
            RackProfileKind::Plate48 => "Cell culture plate 48-well",
            RackProfileKind::Plate96 => "Plate 96",
            RackProfileKind::Plate384 => "Plate 384",
            RackProfileKind::Custom => "Custom",
        }
    }

    pub(super) fn rack_authoring_template_label(template: RackAuthoringTemplate) -> &'static str {
        match template {
            RackAuthoringTemplate::BenchRows => "Bench rows",
            RackAuthoringTemplate::PlateColumns => "Plate columns",
            RackAuthoringTemplate::PlateEdgeAvoidance => "Plate edge avoidance",
        }
    }

    pub(super) fn rack_fill_direction_label(fill_direction: RackFillDirection) -> &'static str {
        match fill_direction {
            RackFillDirection::RowMajor => "Row-major",
            RackFillDirection::ColumnMajor => "Column-major",
        }
    }

    pub(super) fn rack_label_sheet_preset_label(preset: RackLabelSheetPreset) -> &'static str {
        match preset {
            RackLabelSheetPreset::CompactCards => "Compact cards",
            RackLabelSheetPreset::PrintA4 => "Print A4",
            RackLabelSheetPreset::WideCards => "Wide cards",
        }
    }

    pub(super) fn rack_physical_template_label(template: RackPhysicalTemplateKind) -> &'static str {
        match template {
            RackPhysicalTemplateKind::StoragePcrTubeRack => "Storage PCR tube rack",
            RackPhysicalTemplateKind::PipettingPcrTubeRack => "Pipetting PCR tube rack",
            RackPhysicalTemplateKind::CellCulturePlate => "Cell culture plate",
        }
    }

    pub(super) fn rack_carrier_label_preset_label(preset: RackCarrierLabelPreset) -> &'static str {
        match preset {
            RackCarrierLabelPreset::FrontStripAndCards => "Front strip + cards",
            RackCarrierLabelPreset::FrontStripOnly => "Front strip only",
            RackCarrierLabelPreset::ModuleCardsOnly => "Module cards only",
        }
    }

    pub(super) fn rack_coordinate_for_slot(
        profile: &crate::engine::RackProfileSnapshot,
        row: usize,
        column: usize,
    ) -> String {
        GentleEngine::rack_coordinate_from_row_column(profile, row, column).unwrap_or_else(|_| {
            format!(
                "{}{}",
                GentleEngine::rack_row_label_from_index(row),
                column + 1
            )
        })
    }

    pub(super) fn rack_edge_avoidance_coordinates(
        profile: &crate::engine::RackProfileSnapshot,
    ) -> Option<Vec<String>> {
        if profile.rows < 3 || profile.columns < 3 {
            return None;
        }
        let mut out = Vec::new();
        match profile.fill_direction {
            RackFillDirection::RowMajor => {
                for row in 0..profile.rows {
                    for column in 0..profile.columns {
                        if row == 0
                            || row + 1 == profile.rows
                            || column == 0
                            || column + 1 == profile.columns
                        {
                            out.push(Self::rack_coordinate_for_slot(profile, row, column));
                        }
                    }
                }
            }
            RackFillDirection::ColumnMajor => {
                for column in 0..profile.columns {
                    for row in 0..profile.rows {
                        if row == 0
                            || row + 1 == profile.rows
                            || column == 0
                            || column + 1 == profile.columns
                        {
                            out.push(Self::rack_coordinate_for_slot(profile, row, column));
                        }
                    }
                }
            }
        }
        Some(out)
    }

    pub(super) fn infer_rack_authoring_template(
        profile: &crate::engine::RackProfileSnapshot,
    ) -> RackAuthoringTemplate {
        match profile.fill_direction {
            RackFillDirection::RowMajor => RackAuthoringTemplate::BenchRows,
            RackFillDirection::ColumnMajor => {
                if profile.blocked_coordinates.as_slice()
                    == Self::rack_edge_avoidance_coordinates(profile)
                        .unwrap_or_default()
                        .as_slice()
                {
                    RackAuthoringTemplate::PlateEdgeAvoidance
                } else {
                    RackAuthoringTemplate::PlateColumns
                }
            }
        }
    }

    pub(super) fn parse_blocked_coordinate_text(raw: &str) -> Vec<String> {
        raw.split(|ch: char| ch == ',' || ch == ';' || ch.is_whitespace())
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(ToString::to_string)
            .collect()
    }

    pub(super) fn rack_color_for_arrangement(arrangement_id: &str) -> egui::Color32 {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        arrangement_id.hash(&mut hasher);
        let hash = hasher.finish();
        egui::Color32::from_rgb(
            (80 + (hash & 0x3f) as u8).min(220),
            (110 + ((hash >> 8) & 0x5f) as u8).min(225),
            (130 + ((hash >> 16) & 0x5f) as u8).min(235),
        )
    }

    pub(super) fn rack_short_role_label(role_label: &str) -> String {
        let trimmed = role_label.trim();
        if trimmed.eq_ignore_ascii_case("ladder_left") {
            "L".to_string()
        } else if trimmed.eq_ignore_ascii_case("ladder_right") {
            "R".to_string()
        } else if trimmed.eq_ignore_ascii_case("vector") {
            "Vec".to_string()
        } else if trimmed.eq_ignore_ascii_case("product") {
            "Prod".to_string()
        } else if let Some(suffix) = trimmed.strip_prefix("insert_") {
            format!("I{}", suffix)
        } else if let Some(suffix) = trimmed.strip_prefix("lane_") {
            format!("Ln{}", suffix)
        } else {
            trimmed.to_string()
        }
    }

    pub(super) fn arrangement_default_rack_id_from_state(
        &self,
        arrangement_id: &str,
    ) -> Option<String> {
        let engine = self.engine.read().unwrap();
        engine
            .state()
            .container_state
            .arrangements
            .get(arrangement_id.trim())
            .and_then(|arrangement| arrangement.default_rack_id.clone())
            .filter(|rack_id| engine.state().container_state.racks.contains_key(rack_id))
    }

    pub(super) fn ensure_default_rack_for_arrangement_ui(
        &mut self,
        arrangement_id: &str,
    ) -> Option<String> {
        if let Some(rack_id) = self.arrangement_default_rack_id_from_state(arrangement_id) {
            return Some(rack_id);
        }
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::CreateRackFromArrangement {
                arrangement_id: arrangement_id.trim().to_string(),
                rack_id: None,
                name: None,
                profile: None,
            });
        match result {
            Ok(op_result) => {
                self.app_status = op_result.messages.first().cloned().unwrap_or_else(|| {
                    format!(
                        "Created default rack draft for arrangement '{}'",
                        arrangement_id.trim()
                    )
                });
                self.lineage_cache_valid = false;
                self.refresh_lineage_cache_if_needed();
                self.arrangement_default_rack_id_from_state(arrangement_id)
            }
            Err(err) => {
                self.app_status = format!(
                    "Could not create default rack for arrangement '{}': {}",
                    arrangement_id.trim(),
                    err.message
                );
                None
            }
        }
    }

    pub(super) fn open_rack_dialog(&mut self, rack_id: &str) {
        let rack_id = rack_id.trim();
        if rack_id.is_empty() {
            self.app_status = "Rack view requires a non-empty rack id".to_string();
            return;
        }
        self.rack_view_rack_id = rack_id.to_string();
        self.rack_view_status.clear();
        self.rack_view_scroll_offset = Vec2::ZERO;
        self.rack_view_selected_coordinates.clear();
        self.rack_view_selected_arrangement_ids.clear();
        self.rack_view_drag_state = None;
        self.rack_view_hover_target_coordinate = None;
        self.rack_view_recent_drop_ghost = None;
        if let Some(rack) = self
            .engine
            .read()
            .unwrap()
            .state()
            .container_state
            .racks
            .get(rack_id)
            .cloned()
        {
            self.rack_profile_editor_kind = rack.profile.kind;
            self.rack_authoring_template_editor =
                Self::infer_rack_authoring_template(&rack.profile);
            self.rack_fill_direction_editor = rack.profile.fill_direction;
            self.rack_custom_profile_rows = rack.profile.rows.to_string();
            self.rack_custom_profile_columns = rack.profile.columns.to_string();
            self.rack_blocked_coordinates_text = rack.profile.blocked_coordinates.join(", ");
        }
        let was_open = self.show_rack_dialog;
        self.show_rack_dialog = true;
        self.mark_window_open_or_focus(Self::rack_viewport_id(), was_open);
    }

    pub(super) fn open_arrangement_rack_dialog(&mut self, arrangement_id: &str) {
        if let Some(rack_id) = self.ensure_default_rack_for_arrangement_ui(arrangement_id) {
            self.open_rack_dialog(&rack_id);
        }
    }

    pub(super) fn prompt_export_rack_labels_svg_scoped(
        &mut self,
        rack_id: &str,
        arrangement_id: Option<&str>,
        stem: &str,
    ) {
        let rack_id = rack_id.trim();
        if rack_id.is_empty() {
            self.app_status = "Rack label SVG export requires a non-empty rack id".to_string();
            return;
        }
        let arrangement_id = arrangement_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let default_file_name = format!("{stem}.labels.svg");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Rack label SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::ExportRackLabelsSvg {
                rack_id: rack_id.to_string(),
                path: path_text.clone(),
                arrangement_id,
                preset: self.rack_label_sheet_preset,
            });
        match result {
            Ok(op_result) => {
                self.app_status = op_result
                    .messages
                    .first()
                    .cloned()
                    .unwrap_or_else(|| format!("Wrote rack labels SVG to '{path_text}'"));
            }
            Err(err) => {
                self.app_status = format!("Could not export rack labels SVG: {}", err.message);
            }
        }
    }

    pub(super) fn open_arrangement_labels_preview_dialog(&mut self, arrangement_id: &str) {
        let Some(rack_id) = self.ensure_default_rack_for_arrangement_ui(arrangement_id) else {
            return;
        };
        self.open_rack_labels_preview_dialog(&rack_id, Some(arrangement_id));
    }

    pub(super) fn prompt_export_arrangement_carrier_labels_svg(&mut self, arrangement_id: &str) {
        let Some(rack_id) = self.ensure_default_rack_for_arrangement_ui(arrangement_id) else {
            return;
        };
        let stem = Self::sanitize_file_stem(arrangement_id, "rack_carrier_labels");
        let default_file_name = format!("{stem}.carrier.svg");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Rack carrier-label SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::ExportRackCarrierLabelsSvg {
                rack_id: rack_id.clone(),
                path: path_text.clone(),
                arrangement_id: Some(arrangement_id.trim().to_string()),
                template: self.rack_physical_template_kind,
                preset: self.rack_carrier_label_preset,
            });
        match result {
            Ok(op_result) => {
                self.app_status =
                    op_result.messages.first().cloned().unwrap_or_else(|| {
                        format!("Wrote rack carrier-label SVG to '{path_text}'")
                    });
            }
            Err(err) => {
                self.app_status =
                    format!("Could not export rack carrier-label SVG: {}", err.message);
            }
        }
    }

    pub(super) fn prompt_export_rack_fabrication_svg(&mut self, rack_id: &str) {
        let rack_id = rack_id.trim();
        if rack_id.is_empty() {
            self.app_status = "Rack fabrication export requires a non-empty rack id".to_string();
            return;
        }
        let stem = Self::sanitize_file_stem(rack_id, "rack_fabrication");
        let default_file_name = format!("{stem}.fabrication.svg");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Rack fabrication SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::ExportRackFabricationSvg {
                rack_id: rack_id.to_string(),
                path: path_text.clone(),
                template: self.rack_physical_template_kind,
            });
        match result {
            Ok(op_result) => {
                self.app_status = op_result
                    .messages
                    .first()
                    .cloned()
                    .unwrap_or_else(|| format!("Wrote rack fabrication SVG to '{path_text}'"));
                self.rack_view_status = self.app_status.clone();
            }
            Err(err) => {
                self.app_status = format!("Could not export rack fabrication SVG: {}", err.message);
                self.rack_view_status = self.app_status.clone();
            }
        }
    }

    pub(super) fn prompt_export_rack_isometric_svg(&mut self, rack_id: &str) {
        let rack_id = rack_id.trim();
        if rack_id.is_empty() {
            self.app_status = "Rack isometric export requires a non-empty rack id".to_string();
            return;
        }
        let stem = Self::sanitize_file_stem(rack_id, "rack_isometric");
        let default_file_name = format!("{stem}.isometric.svg");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Rack isometric SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::ExportRackIsometricSvg {
                rack_id: rack_id.to_string(),
                path: path_text.clone(),
                template: self.rack_physical_template_kind,
            });
        match result {
            Ok(op_result) => {
                self.app_status = op_result
                    .messages
                    .first()
                    .cloned()
                    .unwrap_or_else(|| format!("Wrote rack isometric SVG to '{path_text}'"));
                self.rack_view_status = self.app_status.clone();
            }
            Err(err) => {
                self.app_status = format!("Could not export rack isometric SVG: {}", err.message);
                self.rack_view_status = self.app_status.clone();
            }
        }
    }

    pub(super) fn prompt_export_rack_openscad(&mut self, rack_id: &str) {
        let rack_id = rack_id.trim();
        if rack_id.is_empty() {
            self.app_status = "Rack OpenSCAD export requires a non-empty rack id".to_string();
            return;
        }
        let stem = Self::sanitize_file_stem(rack_id, "rack");
        let default_file_name = format!("{stem}.scad");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("OpenSCAD", &["scad"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Rack OpenSCAD export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::ExportRackOpenScad {
                rack_id: rack_id.to_string(),
                path: path_text.clone(),
                template: self.rack_physical_template_kind,
            });
        match result {
            Ok(op_result) => {
                self.app_status = op_result
                    .messages
                    .first()
                    .cloned()
                    .unwrap_or_else(|| format!("Wrote rack OpenSCAD to '{path_text}'"));
                self.rack_view_status = self.app_status.clone();
            }
            Err(err) => {
                self.app_status = format!("Could not export rack OpenSCAD: {}", err.message);
                self.rack_view_status = self.app_status.clone();
            }
        }
    }

    pub(super) fn prompt_export_rack_carrier_labels_svg(&mut self, rack_id: &str) {
        let rack_id = rack_id.trim();
        if rack_id.is_empty() {
            self.app_status = "Rack carrier-label export requires a non-empty rack id".to_string();
            return;
        }
        let stem = Self::sanitize_file_stem(rack_id, "rack_carrier_labels");
        let default_file_name = format!("{stem}.carrier.svg");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Rack carrier-label SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::ExportRackCarrierLabelsSvg {
                rack_id: rack_id.to_string(),
                path: path_text.clone(),
                arrangement_id: None,
                template: self.rack_physical_template_kind,
                preset: self.rack_carrier_label_preset,
            });
        match result {
            Ok(op_result) => {
                self.app_status =
                    op_result.messages.first().cloned().unwrap_or_else(|| {
                        format!("Wrote rack carrier-label SVG to '{path_text}'")
                    });
                self.rack_view_status = self.app_status.clone();
            }
            Err(err) => {
                self.app_status =
                    format!("Could not export rack carrier-label SVG: {}", err.message);
                self.rack_view_status = self.app_status.clone();
            }
        }
    }

    pub(super) fn prompt_export_rack_simulation_json(&mut self, rack_id: &str) {
        let rack_id = rack_id.trim();
        if rack_id.is_empty() {
            self.app_status = "Rack simulation export requires a non-empty rack id".to_string();
            return;
        }
        let stem = Self::sanitize_file_stem(rack_id, "rack_simulation");
        let default_file_name = format!("{stem}.simulation.json");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("JSON", &["json"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Rack simulation JSON export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::ExportRackSimulationJson {
                rack_id: rack_id.to_string(),
                path: path_text.clone(),
                template: self.rack_physical_template_kind,
            });
        match result {
            Ok(op_result) => {
                self.app_status = op_result
                    .messages
                    .first()
                    .cloned()
                    .unwrap_or_else(|| format!("Wrote rack simulation JSON to '{path_text}'"));
                self.rack_view_status = self.app_status.clone();
            }
            Err(err) => {
                self.app_status = format!("Could not export rack simulation JSON: {}", err.message);
                self.rack_view_status = self.app_status.clone();
            }
        }
    }

    pub(super) fn current_rack_snapshot(&self) -> Option<Rack> {
        self.rack_snapshot_by_id(&self.rack_view_rack_id)
    }

    pub(super) fn rack_snapshot_by_id(&self, rack_id: &str) -> Option<Rack> {
        let rack_id = rack_id.trim();
        if rack_id.is_empty() {
            return None;
        }
        self.engine
            .read()
            .unwrap()
            .state()
            .container_state
            .racks
            .get(rack_id)
            .cloned()
    }

    pub(super) fn rack_sorted_entries(
        rack: &Rack,
    ) -> Vec<(usize, String, crate::engine::RackPlacementEntry)> {
        let mut out = rack
            .placements
            .iter()
            .filter_map(|entry| {
                GentleEngine::rack_index_from_coordinate(&rack.profile, &entry.coordinate)
                    .ok()
                    .map(|index| (index, entry.coordinate.clone(), entry.clone()))
            })
            .collect::<Vec<_>>();
        out.sort_by(|a, b| a.0.cmp(&b.0).then(a.2.order_index.cmp(&b.2.order_index)));
        out
    }

    pub(super) fn rack_first_coordinate_for_arrangement(
        sorted_entries: &[(usize, String, crate::engine::RackPlacementEntry)],
        arrangement_id: &str,
    ) -> Option<String> {
        sorted_entries
            .iter()
            .filter(|(_, _, entry)| entry.arrangement_id == arrangement_id)
            .min_by_key(|(index, _, _)| *index)
            .map(|(_, coordinate, _)| coordinate.clone())
    }

    pub(super) fn rack_selected_coordinates_in_order(
        sorted_entries: &[(usize, String, crate::engine::RackPlacementEntry)],
        selected_coordinates: &BTreeSet<String>,
    ) -> Vec<String> {
        sorted_entries
            .iter()
            .filter_map(|(_, coordinate, _)| {
                if selected_coordinates.contains(coordinate) {
                    Some(coordinate.clone())
                } else {
                    None
                }
            })
            .collect::<Vec<_>>()
    }

    pub(super) fn rack_selected_samples_arrangement_id(
        sorted_entries: &[(usize, String, crate::engine::RackPlacementEntry)],
        selected_coordinates: &BTreeSet<String>,
    ) -> Option<String> {
        let mut arrangement_id: Option<String> = None;
        for (_, coordinate, entry) in sorted_entries {
            if selected_coordinates.contains(coordinate) {
                match arrangement_id.as_deref() {
                    Some(existing) if existing != entry.arrangement_id => return None,
                    Some(_) => {}
                    None => arrangement_id = Some(entry.arrangement_id.clone()),
                }
            }
        }
        arrangement_id
    }

    pub(super) fn rack_selected_arrangement_ids_in_order(
        sorted_entries: &[(usize, String, crate::engine::RackPlacementEntry)],
        selected_arrangement_ids: &BTreeSet<String>,
    ) -> Vec<String> {
        let mut seen = HashSet::new();
        let mut out = Vec::new();
        for (_, _, entry) in sorted_entries {
            if selected_arrangement_ids.contains(&entry.arrangement_id)
                && seen.insert(entry.arrangement_id.clone())
            {
                out.push(entry.arrangement_id.clone());
            }
        }
        out
    }

    pub(super) fn rack_drag_status_text(drag: &RackDragState, drop_target: Option<&str>) -> String {
        match drag {
            RackDragState::Sample {
                from_coordinate,
                arrangement_id,
                role_label,
            } => match drop_target {
                Some(target) => format!(
                    "Dragging sample '{}' from {} in arrangement '{}' to {}. Release to shift neighboring slots.",
                    role_label, from_coordinate, arrangement_id, target
                ),
                None => format!(
                    "Dragging sample '{}' from {} in arrangement '{}'. Hover a rack slot and release to shift neighboring slots.",
                    role_label, from_coordinate, arrangement_id
                ),
            },
            RackDragState::Samples {
                arrangement_id,
                from_coordinates,
            } => {
                let label = if from_coordinates.len() <= 3 {
                    from_coordinates.join(", ")
                } else {
                    format!(
                        "{} (+{} more)",
                        from_coordinates[..3].join(", "),
                        from_coordinates.len() - 3
                    )
                };
                match drop_target {
                    Some(target) => format!(
                        "Dragging {} samples [{}] in arrangement '{}' to {}. Release to move the selected samples together within that arrangement block.",
                        from_coordinates.len(),
                        label,
                        arrangement_id,
                        target
                    ),
                    None => format!(
                        "Dragging {} samples [{}] in arrangement '{}'. Hover a rack slot and release to move the selected samples together.",
                        from_coordinates.len(),
                        label,
                        arrangement_id
                    ),
                }
            }
            RackDragState::ArrangementBlock {
                arrangement_id,
                from_coordinate,
            } => match drop_target {
                Some(target) => format!(
                    "Dragging arrangement block '{}' from {} to {}. Release to move the whole block and shift later occupied slots.",
                    arrangement_id, from_coordinate, target
                ),
                None => format!(
                    "Dragging arrangement block '{}' from {}. Hover a rack slot and release to move the whole block.",
                    arrangement_id, from_coordinate
                ),
            },
            RackDragState::ArrangementBlocks {
                arrangement_ids,
                from_coordinate,
            } => {
                let label = if arrangement_ids.len() <= 3 {
                    arrangement_ids.join(", ")
                } else {
                    format!(
                        "{} (+{} more)",
                        arrangement_ids[..3].join(", "),
                        arrangement_ids.len() - 3
                    )
                };
                match drop_target {
                    Some(target) => format!(
                        "Dragging {} arrangement blocks [{}] from {} to {}. Release to move the selected blocks together and shift later occupied slots.",
                        arrangement_ids.len(),
                        label,
                        from_coordinate,
                        target
                    ),
                    None => format!(
                        "Dragging {} arrangement blocks [{}] from {}. Hover a rack slot and release to move the selected blocks together.",
                        arrangement_ids.len(),
                        label,
                        from_coordinate
                    ),
                }
            }
        }
    }

    pub(super) fn rack_help_strip_items() -> &'static [(&'static str, &'static str)] {
        &[
            (
                "Samples",
                "Drag one slot, or multi-select same-arrangement samples to move them together.",
            ),
            (
                "Blocks",
                "Drag arrangement chips, or multi-select chips to move whole experiment blocks.",
            ),
            (
                "Preview",
                "Pulsing ghosts show the predicted drop order; changed slots fade briefly after drop.",
            ),
            (
                "Keys",
                "Cancels the active drag, and click-select plus click-target still works.",
            ),
        ]
    }

    pub(super) fn rack_help_strip_keycaps(title: &str) -> &'static [&'static str] {
        match title {
            "Samples" | "Blocks" => &["Cmd", "Ctrl"],
            "Keys" => &["Esc"],
            _ => &[],
        }
    }

    pub(super) fn render_rack_help_keycap(ui: &mut Ui, label: &str) {
        egui::Frame::group(ui.style())
            .fill(ui.visuals().faint_bg_color)
            .stroke(egui::Stroke::new(
                1.0,
                ui.visuals().widgets.noninteractive.bg_stroke.color,
            ))
            .inner_margin(egui::Margin::symmetric(5, 2))
            .show(ui, |ui| {
                ui.small(egui::RichText::new(label).strong());
            });
    }

    pub(super) fn rack_help_toggle_label(collapsed: bool) -> &'static str {
        if collapsed { "Show help" } else { "Hide help" }
    }

    pub(super) fn rack_help_pin_label(pinned: bool) -> &'static str {
        if pinned { "Unpin" } else { "Pin open" }
    }

    pub(super) fn record_successful_rack_move_and_maybe_autocollapse(&mut self) {
        if self.rack_help_strip_pinned_open || self.rack_help_strip_auto_minimized {
            return;
        }
        self.rack_help_strip_successful_move_count =
            self.rack_help_strip_successful_move_count.saturating_add(1);
        if self.rack_help_strip_successful_move_count >= RACK_HELP_AUTO_MINIMIZE_MOVE_THRESHOLD {
            self.rack_help_strip_auto_minimized = true;
            self.rack_help_strip_collapsed = true;
            self.persist_rack_workspace_to_state();
            self.rack_view_status = format!(
                "Rack help auto-minimized after {} successful moves. Use Show help to reopen or Pin open to keep it visible.",
                RACK_HELP_AUTO_MINIMIZE_MOVE_THRESHOLD
            );
        } else {
            self.persist_rack_workspace_to_state();
        }
    }

    pub(super) fn render_rack_help_strip(&mut self, ui: &mut Ui, arrangement_ids: &[String]) {
        egui::Frame::group(ui.style())
            .inner_margin(egui::Margin::symmetric(8, 6))
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.small(egui::RichText::new("Rack help").strong());
                    let toggle = Self::rack_help_toggle_label(self.rack_help_strip_collapsed);
                    if ui
                        .small_button(toggle)
                        .on_hover_text("Collapse or reopen the rack interaction help strip")
                        .clicked()
                    {
                        self.rack_help_strip_collapsed = !self.rack_help_strip_collapsed;
                        self.persist_rack_workspace_to_state();
                    }
                    let pin_label = Self::rack_help_pin_label(self.rack_help_strip_pinned_open);
                    if ui
                        .small_button(pin_label)
                        .on_hover_text(
                            "Keep the rack help strip open instead of auto-minimizing after a few successful rack moves",
                        )
                        .clicked()
                    {
                        self.rack_help_strip_pinned_open = !self.rack_help_strip_pinned_open;
                        if self.rack_help_strip_pinned_open {
                            self.rack_help_strip_collapsed = false;
                        }
                        self.persist_rack_workspace_to_state();
                    }
                    if !self.rack_help_strip_pinned_open && !self.rack_help_strip_auto_minimized {
                        ui.small(format!(
                            "auto-hide after {} moves",
                            RACK_HELP_AUTO_MINIMIZE_MOVE_THRESHOLD
                        ));
                    } else if self.rack_help_strip_pinned_open {
                        ui.small("pinned open");
                    }
                });
                if self.rack_help_strip_collapsed {
                    return;
                }
                ui.add_space(4.0);
                ui.horizontal_wrapped(|ui| {
                    if !arrangement_ids.is_empty() {
                        ui.small("Arrangement colors");
                        for arrangement_id in arrangement_ids.iter().take(4) {
                            let color = Self::rack_color_for_arrangement(arrangement_id);
                            egui::Frame::group(ui.style())
                                .fill(color.gamma_multiply(0.14))
                                .stroke(egui::Stroke::new(1.0, color.gamma_multiply(0.85)))
                                .inner_margin(egui::Margin::symmetric(6, 3))
                                .show(ui, |ui| {
                                    ui.small(
                                        egui::RichText::new(arrangement_id.clone())
                                            .color(color)
                                            .strong(),
                                    );
                                });
                        }
                        if arrangement_ids.len() > 4 {
                            ui.small(format!("+{} more", arrangement_ids.len() - 4));
                        }
                    }
                });
                ui.add_space(4.0);
                ui.horizontal_wrapped(|ui| {
                    for (idx, (title, body)) in Self::rack_help_strip_items().iter().enumerate() {
                        let chip_color = arrangement_ids
                            .get(idx % arrangement_ids.len().max(1))
                            .map(|arrangement_id| Self::rack_color_for_arrangement(arrangement_id))
                            .unwrap_or(ui.visuals().strong_text_color());
                        ui.group(|ui| {
                            ui.set_max_width(220.0);
                            ui.horizontal(|ui| {
                                ui.small(egui::RichText::new("■").color(chip_color));
                                ui.small(egui::RichText::new(*title).strong());
                            });
                            let keycaps = Self::rack_help_strip_keycaps(title);
                            if !keycaps.is_empty() {
                                ui.horizontal_wrapped(|ui| {
                                    for keycap in keycaps {
                                        Self::render_rack_help_keycap(ui, keycap);
                                    }
                                });
                            }
                            ui.small(*body);
                        });
                    }
                });
            });
    }

    pub(super) fn rack_available_coordinates(
        profile: &crate::engine::RackProfileSnapshot,
    ) -> Vec<String> {
        let blocked = profile
            .blocked_coordinates
            .iter()
            .cloned()
            .collect::<HashSet<_>>();
        let mut coordinates = Vec::new();
        match profile.fill_direction {
            RackFillDirection::RowMajor => {
                for row in 0..profile.rows {
                    for column in 0..profile.columns {
                        let coordinate = Self::rack_coordinate_for_slot(profile, row, column);
                        if !blocked.contains(&coordinate) {
                            coordinates.push(coordinate);
                        }
                    }
                }
            }
            RackFillDirection::ColumnMajor => {
                for column in 0..profile.columns {
                    for row in 0..profile.rows {
                        let coordinate = Self::rack_coordinate_for_slot(profile, row, column);
                        if !blocked.contains(&coordinate) {
                            coordinates.push(coordinate);
                        }
                    }
                }
            }
        }
        coordinates
    }

    pub(super) fn rack_block_insertion_index(
        ordered: &[(usize, crate::engine::RackPlacementEntry)],
        target_index: usize,
    ) -> usize {
        let mut idx = 0usize;
        while idx < ordered.len() {
            let block_start_index = ordered[idx].0;
            if target_index <= block_start_index {
                return idx;
            }
            let arrangement_id = ordered[idx].1.arrangement_id.clone();
            while idx < ordered.len() && ordered[idx].1.arrangement_id == arrangement_id {
                idx += 1;
            }
        }
        ordered.len()
    }

    pub(super) fn rack_block_insertion_index_after_target(
        ordered: &[(usize, crate::engine::RackPlacementEntry)],
        target_index: usize,
    ) -> usize {
        let mut idx = 0usize;
        while idx < ordered.len() {
            let block_start_index = ordered[idx].0;
            if target_index < block_start_index {
                return idx;
            }
            let arrangement_id = ordered[idx].1.arrangement_id.clone();
            let mut block_end_index = block_start_index;
            while idx < ordered.len() && ordered[idx].1.arrangement_id == arrangement_id {
                block_end_index = ordered[idx].0;
                idx += 1;
            }
            if target_index <= block_end_index {
                return idx;
            }
        }
        ordered.len()
    }

    pub(super) fn rack_group_insertion_index_within_block(
        local_to: usize,
        remaining_len: usize,
    ) -> usize {
        local_to.min(remaining_len)
    }

    pub(super) fn rack_changed_coordinates(before: &Rack, after: &Rack) -> BTreeSet<String> {
        let before_map = before
            .placements
            .iter()
            .map(|entry| (entry.coordinate.clone(), entry.clone()))
            .collect::<BTreeMap<_, _>>();
        let after_map = after
            .placements
            .iter()
            .map(|entry| (entry.coordinate.clone(), entry.clone()))
            .collect::<BTreeMap<_, _>>();
        before_map
            .keys()
            .chain(after_map.keys())
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .filter(|coordinate| before_map.get(coordinate) != after_map.get(coordinate))
            .collect::<BTreeSet<_>>()
    }

    pub(super) fn rack_preview_entries_for_drop(
        rack: &Rack,
        sorted_entries: &[(usize, String, crate::engine::RackPlacementEntry)],
        drag: &RackDragState,
        drop_target: &str,
    ) -> Option<Vec<crate::engine::RackPlacementEntry>> {
        let available_coordinates = Self::rack_available_coordinates(&rack.profile);
        let mut ordered = sorted_entries
            .iter()
            .map(|(index, _, entry)| (*index, entry.clone()))
            .collect::<Vec<_>>();
        let from_coordinate = match drag {
            RackDragState::Sample {
                from_coordinate, ..
            } => from_coordinate.as_str(),
            RackDragState::Samples {
                from_coordinates, ..
            } => from_coordinates.first()?.as_str(),
            RackDragState::ArrangementBlock {
                from_coordinate, ..
            }
            | RackDragState::ArrangementBlocks {
                from_coordinate, ..
            } => from_coordinate.as_str(),
        };
        let from_index =
            GentleEngine::rack_index_from_coordinate(&rack.profile, from_coordinate).ok()?;
        let from_pos = ordered.iter().position(|(index, _)| *index == from_index)?;
        match drag {
            RackDragState::ArrangementBlocks {
                arrangement_ids, ..
            } => {
                let to_index =
                    GentleEngine::rack_index_from_coordinate(&rack.profile, drop_target).ok()?;
                let selected = arrangement_ids.iter().cloned().collect::<HashSet<_>>();
                let selected_entries = ordered
                    .iter()
                    .filter(|(_, entry)| selected.contains(&entry.arrangement_id))
                    .map(|(_, entry)| entry.clone())
                    .collect::<Vec<_>>();
                ordered.retain(|(_, entry)| !selected.contains(&entry.arrangement_id));
                if ordered.len() + selected_entries.len() > available_coordinates.len() {
                    return None;
                }
                let insertion_index = Self::rack_block_insertion_index(&ordered, to_index);
                for (offset, entry) in selected_entries.into_iter().enumerate() {
                    ordered.insert(insertion_index + offset, (insertion_index + offset, entry));
                }
            }
            RackDragState::ArrangementBlock { arrangement_id, .. } => {
                let to_index =
                    GentleEngine::rack_index_from_coordinate(&rack.profile, drop_target).ok()?;
                let block_entries = ordered
                    .iter()
                    .filter(|(_, entry)| entry.arrangement_id == *arrangement_id)
                    .map(|(_, entry)| entry.clone())
                    .collect::<Vec<_>>();
                ordered.retain(|(_, entry)| entry.arrangement_id != *arrangement_id);
                if ordered.len() + block_entries.len() > available_coordinates.len() {
                    return None;
                }
                let insertion_index =
                    Self::rack_block_insertion_index_after_target(&ordered, to_index);
                for (offset, entry) in block_entries.into_iter().enumerate() {
                    ordered.insert(insertion_index + offset, (insertion_index + offset, entry));
                }
            }
            RackDragState::Samples {
                arrangement_id,
                from_coordinates,
            } => {
                let target_pos = sorted_entries
                    .iter()
                    .position(|(_, coordinate, _)| coordinate == drop_target)?;
                let block_positions = ordered
                    .iter()
                    .enumerate()
                    .filter_map(|(pos, (_, entry))| {
                        if entry.arrangement_id == *arrangement_id {
                            Some(pos)
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();
                let block_start = *block_positions.first().unwrap_or(&from_pos);
                let block_end = *block_positions.last().unwrap_or(&from_pos);
                if target_pos < block_start || target_pos > block_end {
                    return None;
                }
                let mut block_entries = ordered[block_start..=block_end]
                    .iter()
                    .map(|(_, entry)| entry.clone())
                    .collect::<Vec<_>>();
                let selected_coordinates =
                    from_coordinates.iter().cloned().collect::<BTreeSet<_>>();
                let selected_local_positions = ordered[block_start..=block_end]
                    .iter()
                    .enumerate()
                    .filter_map(|(pos, (_, entry))| {
                        if selected_coordinates.contains(&entry.coordinate) {
                            Some(pos)
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();
                let selected_local_set = selected_local_positions
                    .iter()
                    .copied()
                    .collect::<BTreeSet<_>>();
                let selected_entries = block_entries
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, entry)| {
                        if selected_local_set.contains(&idx) {
                            Some(entry.clone())
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();
                let remaining_entries = block_entries
                    .drain(..)
                    .enumerate()
                    .filter_map(|(idx, entry)| {
                        if selected_local_set.contains(&idx) {
                            None
                        } else {
                            Some(entry)
                        }
                    })
                    .collect::<Vec<_>>();
                let local_to = target_pos - block_start;
                let insertion_index = Self::rack_group_insertion_index_within_block(
                    local_to,
                    remaining_entries.len(),
                );
                let mut reordered_block = Vec::with_capacity(block_end - block_start + 1);
                reordered_block.extend(remaining_entries[..insertion_index].iter().cloned());
                reordered_block.extend(selected_entries);
                reordered_block.extend(remaining_entries[insertion_index..].iter().cloned());
                for (offset, entry) in reordered_block.into_iter().enumerate() {
                    ordered[block_start + offset].1 = entry;
                }
            }
            RackDragState::Sample { arrangement_id, .. } => {
                let target_pos = sorted_entries
                    .iter()
                    .position(|(_, coordinate, _)| coordinate == drop_target)?;
                let block_positions = ordered
                    .iter()
                    .enumerate()
                    .filter_map(|(pos, (_, entry))| {
                        if entry.arrangement_id == *arrangement_id {
                            Some(pos)
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();
                let block_start = *block_positions.first().unwrap_or(&from_pos);
                let block_end = *block_positions.last().unwrap_or(&from_pos);
                if target_pos < block_start || target_pos > block_end {
                    return None;
                }
                let mut block_entries = ordered[block_start..=block_end]
                    .iter()
                    .map(|(_, entry)| entry.clone())
                    .collect::<Vec<_>>();
                let local_from = from_pos - block_start;
                let local_to = target_pos - block_start;
                let moved = block_entries.remove(local_from);
                block_entries.insert(local_to, moved);
                for (offset, entry) in block_entries.into_iter().enumerate() {
                    ordered[block_start + offset].1 = entry;
                }
            }
        }
        let mut reflowed = ordered
            .into_iter()
            .map(|(_, entry)| entry)
            .collect::<Vec<_>>();
        if reflowed.len() > available_coordinates.len() {
            return None;
        }
        let mut order_by_arrangement: HashMap<String, usize> = HashMap::new();
        for (entry, coordinate) in reflowed.iter_mut().zip(available_coordinates) {
            entry.coordinate = coordinate;
            let next = order_by_arrangement
                .entry(entry.arrangement_id.clone())
                .and_modify(|index| *index += 1)
                .or_insert(0usize);
            entry.order_index = *next;
        }
        Some(reflowed)
    }

    pub(super) fn rack_autoscroll_offset_for_pointer(
        current_offset: Vec2,
        pointer: Pos2,
        viewport_rect: egui::Rect,
        max_scroll_x: f32,
        max_scroll_y: f32,
    ) -> Vec2 {
        const THRESHOLD_PX: f32 = 36.0;
        const MAX_STEP_PX: f32 = 24.0;
        let mut next = current_offset;
        let left_zone = viewport_rect.left() + THRESHOLD_PX;
        let right_zone = viewport_rect.right() - THRESHOLD_PX;
        let top_zone = viewport_rect.top() + THRESHOLD_PX;
        let bottom_zone = viewport_rect.bottom() - THRESHOLD_PX;
        if pointer.x < left_zone {
            let factor = ((left_zone - pointer.x) / THRESHOLD_PX).clamp(0.0, 1.0);
            next.x = (next.x - factor * MAX_STEP_PX).max(0.0);
        } else if pointer.x > right_zone {
            let factor = ((pointer.x - right_zone) / THRESHOLD_PX).clamp(0.0, 1.0);
            next.x = (next.x + factor * MAX_STEP_PX).min(max_scroll_x);
        }
        if pointer.y < top_zone {
            let factor = ((top_zone - pointer.y) / THRESHOLD_PX).clamp(0.0, 1.0);
            next.y = (next.y - factor * MAX_STEP_PX).max(0.0);
        } else if pointer.y > bottom_zone {
            let factor = ((pointer.y - bottom_zone) / THRESHOLD_PX).clamp(0.0, 1.0);
            next.y = (next.y + factor * MAX_STEP_PX).min(max_scroll_y);
        }
        next
    }

    pub(super) fn rack_ghost_preview_map(
        rack: &Rack,
        sorted_entries: &[(usize, String, crate::engine::RackPlacementEntry)],
        drag: &RackDragState,
        drop_target: &str,
    ) -> Option<BTreeMap<String, RackGhostPreviewCell>> {
        let preview_entries =
            Self::rack_preview_entries_for_drop(rack, sorted_entries, drag, drop_target)?;
        let current_entries = sorted_entries
            .iter()
            .map(|(_, coordinate, entry)| (coordinate.clone(), entry.clone()))
            .collect::<BTreeMap<_, _>>();
        let preview_map = preview_entries
            .into_iter()
            .map(|entry| (entry.coordinate.clone(), entry))
            .collect::<BTreeMap<_, _>>();
        let mut out = BTreeMap::new();
        let all_coordinates = current_entries
            .keys()
            .chain(preview_map.keys())
            .cloned()
            .collect::<BTreeSet<_>>();
        for coordinate in all_coordinates {
            let current_entry = current_entries.get(&coordinate);
            let preview_entry = preview_map.get(&coordinate);
            if current_entry != preview_entry {
                out.insert(
                    coordinate,
                    RackGhostPreviewCell {
                        predicted_entry: preview_entry.cloned(),
                    },
                );
            }
        }
        Some(out)
    }

    pub(super) fn rack_drag_preview_animation_alpha(now_seconds: f64) -> f32 {
        let oscillation = (((now_seconds * 5.0).sin() as f32) + 1.0) * 0.5;
        0.08 + oscillation * 0.08
    }

    pub(super) fn record_rack_drop_ghost(&mut self, rack_id: &str, before_rack: Option<Rack>) {
        let Some(before_rack) = before_rack else {
            self.rack_view_recent_drop_ghost = None;
            return;
        };
        let Some(after_rack) = self.rack_snapshot_by_id(rack_id) else {
            self.rack_view_recent_drop_ghost = None;
            return;
        };
        let changed_coordinates = Self::rack_changed_coordinates(&before_rack, &after_rack);
        if changed_coordinates.is_empty() {
            self.rack_view_recent_drop_ghost = None;
        } else {
            self.rack_view_recent_drop_ghost = Some(RackDropGhostState {
                rack_id: rack_id.trim().to_string(),
                changed_coordinates,
                started_at: Instant::now(),
            });
        }
    }

    pub(super) fn rack_entry_display(entry: &crate::engine::RackPlacementEntry) -> String {
        let role_text = Self::rack_short_role_label(&entry.role_label);
        match entry.occupant.as_ref() {
            Some(RackOccupant::Container { container_id }) => {
                format!("{role_text}\n{container_id}")
            }
            Some(RackOccupant::LadderReference { ladder_name }) => {
                format!("{role_text}\n{ladder_name}")
            }
            None => role_text,
        }
    }

    pub(super) fn apply_rack_move(
        &mut self,
        rack_id: &str,
        from_coordinate: &str,
        to_coordinate: &str,
        move_block: bool,
    ) {
        if from_coordinate.trim().is_empty() || to_coordinate.trim().is_empty() {
            self.rack_view_status =
                "Rack move requires both a source and target coordinate.".to_string();
            return;
        }
        if from_coordinate.trim() == to_coordinate.trim() {
            self.rack_view_status =
                "Rack move canceled because the source and target are the same slot.".to_string();
            return;
        }
        let before_rack = self.rack_snapshot_by_id(rack_id);
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::MoveRackPlacement {
                rack_id: rack_id.trim().to_string(),
                from_coordinate: from_coordinate.trim().to_string(),
                to_coordinate: to_coordinate.trim().to_string(),
                move_block,
            });
        match result {
            Ok(op_result) => {
                self.rack_view_status = op_result.messages.first().cloned().unwrap_or_else(|| {
                    if move_block {
                        "Moved arrangement block".to_string()
                    } else {
                        "Moved sample on rack".to_string()
                    }
                });
                self.lineage_cache_valid = false;
                self.refresh_lineage_cache_if_needed();
                self.record_rack_drop_ghost(rack_id, before_rack);
                self.record_successful_rack_move_and_maybe_autocollapse();
            }
            Err(err) => {
                self.rack_view_status = if move_block {
                    format!("Could not move arrangement block: {}", err.message)
                } else {
                    format!("Could not move sample: {}", err.message)
                };
            }
        }
    }

    pub(super) fn apply_rack_move_blocks(
        &mut self,
        rack_id: &str,
        arrangement_ids: &[String],
        to_coordinate: &str,
    ) {
        if arrangement_ids.is_empty() || to_coordinate.trim().is_empty() {
            self.rack_view_status =
                "Rack block move requires selected arrangement blocks and a target coordinate."
                    .to_string();
            return;
        }
        let before_rack = self.rack_snapshot_by_id(rack_id);
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::MoveRackArrangementBlocks {
                rack_id: rack_id.trim().to_string(),
                arrangement_ids: arrangement_ids.to_vec(),
                to_coordinate: to_coordinate.trim().to_string(),
            });
        match result {
            Ok(op_result) => {
                self.rack_view_status = op_result.messages.first().cloned().unwrap_or_else(|| {
                    format!("Moved {} arrangement blocks on rack", arrangement_ids.len())
                });
                self.lineage_cache_valid = false;
                self.refresh_lineage_cache_if_needed();
                self.record_rack_drop_ghost(rack_id, before_rack);
                self.record_successful_rack_move_and_maybe_autocollapse();
            }
            Err(err) => {
                self.rack_view_status =
                    format!("Could not move arrangement blocks: {}", err.message);
            }
        }
    }

    pub(super) fn apply_rack_move_samples(
        &mut self,
        rack_id: &str,
        from_coordinates: &[String],
        to_coordinate: &str,
    ) {
        if from_coordinates.is_empty() || to_coordinate.trim().is_empty() {
            self.rack_view_status =
                "Rack sample-group move requires selected samples and a target coordinate."
                    .to_string();
            return;
        }
        let before_rack = self.rack_snapshot_by_id(rack_id);
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::MoveRackSamples {
                rack_id: rack_id.trim().to_string(),
                from_coordinates: from_coordinates.to_vec(),
                to_coordinate: to_coordinate.trim().to_string(),
            });
        match result {
            Ok(op_result) => {
                self.rack_view_status =
                    op_result.messages.first().cloned().unwrap_or_else(|| {
                        format!("Moved {} samples on rack", from_coordinates.len())
                    });
                self.lineage_cache_valid = false;
                self.refresh_lineage_cache_if_needed();
                self.record_rack_drop_ghost(rack_id, before_rack);
                self.record_successful_rack_move_and_maybe_autocollapse();
            }
            Err(err) => {
                self.rack_view_status = format!("Could not move selected samples: {}", err.message);
            }
        }
    }

    pub(super) fn render_rack_contents(&mut self, ui: &mut Ui) -> bool {
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text("Rack");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        let Some(rack) = self.current_rack_snapshot() else {
            ui.small("No rack is currently selected.");
            return close_requested;
        };
        ui.label("Physical rack/plate placement linked to one or more arrangements. Drag one occupied sample slot to another slot to shift neighbors within that arrangement block, or drag an arrangement chip to move the whole block. Command/Ctrl-click arrangement chips to select multiple blocks, then drag one selected chip or click a target slot. Click-select and click-target still works as a fallback.");
        ui.horizontal(|ui| {
            ui.label("Rack");
            ui.monospace(format!("{} ({})", rack.name, rack.rack_id));
            ui.label("Profile");
            let previous_profile_kind = self.rack_profile_editor_kind;
            egui::ComboBox::from_id_salt("rack_profile_combo")
                .selected_text(Self::rack_profile_label(self.rack_profile_editor_kind))
                .show_ui(ui, |ui| {
                    for profile in [
                        RackProfileKind::SmallTube4x6,
                        RackProfileKind::Plate6,
                        RackProfileKind::Plate12,
                        RackProfileKind::Plate24,
                        RackProfileKind::Plate48,
                        RackProfileKind::Plate96,
                        RackProfileKind::Plate384,
                        RackProfileKind::Custom,
                    ] {
                        ui.selectable_value(
                            &mut self.rack_profile_editor_kind,
                            profile,
                            Self::rack_profile_label(profile),
                        );
                    }
                });
            ui.label("Fill");
            let previous_fill_direction = self.rack_fill_direction_editor;
            egui::ComboBox::from_id_salt("rack_fill_direction_combo")
                .selected_text(Self::rack_fill_direction_label(
                    self.rack_fill_direction_editor,
                ))
                .show_ui(ui, |ui| {
                    for fill_direction in
                        [RackFillDirection::RowMajor, RackFillDirection::ColumnMajor]
                    {
                        ui.selectable_value(
                            &mut self.rack_fill_direction_editor,
                            fill_direction,
                            Self::rack_fill_direction_label(fill_direction),
                        );
                    }
                });
            if self.rack_profile_editor_kind != previous_profile_kind
                && self.rack_profile_editor_kind == RackProfileKind::Custom
            {
                self.rack_custom_profile_rows = rack.profile.rows.to_string();
                self.rack_custom_profile_columns = rack.profile.columns.to_string();
                self.rack_view_status =
                    "Custom rack profile selected. Edit rows/columns and click Apply Custom."
                        .to_string();
            }
            if self.rack_profile_editor_kind != rack.profile.kind
                && self.rack_profile_editor_kind != RackProfileKind::Custom
            {
                let result = self
                    .engine
                    .write()
                    .unwrap()
                    .apply(Operation::SetRackProfile {
                        rack_id: rack.rack_id.clone(),
                        profile: self.rack_profile_editor_kind,
                    });
                match result {
                    Ok(op_result) => {
                        self.rack_view_status = op_result
                            .messages
                            .first()
                            .cloned()
                            .unwrap_or_else(|| "Updated rack profile".to_string());
                        self.lineage_cache_valid = false;
                        self.refresh_lineage_cache_if_needed();
                        if let Some(updated_rack) = self.current_rack_snapshot() {
                            self.rack_custom_profile_rows = updated_rack.profile.rows.to_string();
                            self.rack_custom_profile_columns =
                                updated_rack.profile.columns.to_string();
                            self.rack_profile_editor_kind = updated_rack.profile.kind;
                            self.rack_authoring_template_editor =
                                Self::infer_rack_authoring_template(&updated_rack.profile);
                            self.rack_fill_direction_editor = updated_rack.profile.fill_direction;
                            self.rack_blocked_coordinates_text =
                                updated_rack.profile.blocked_coordinates.join(", ");
                        }
                    }
                    Err(err) => {
                        self.rack_view_status =
                            format!("Could not change rack profile: {}", err.message);
                        self.rack_profile_editor_kind = rack.profile.kind;
                    }
                }
            }
            if self.rack_fill_direction_editor != previous_fill_direction {
                let result = self
                    .engine
                    .write()
                    .unwrap()
                    .apply(Operation::SetRackFillDirection {
                        rack_id: rack.rack_id.clone(),
                        fill_direction: self.rack_fill_direction_editor,
                    });
                match result {
                    Ok(op_result) => {
                        self.rack_view_status = op_result
                            .messages
                            .first()
                            .cloned()
                            .unwrap_or_else(|| "Updated rack fill direction".to_string());
                        self.lineage_cache_valid = false;
                        self.refresh_lineage_cache_if_needed();
                        if let Some(updated_rack) = self.current_rack_snapshot() {
                            self.rack_authoring_template_editor =
                                Self::infer_rack_authoring_template(&updated_rack.profile);
                            self.rack_fill_direction_editor = updated_rack.profile.fill_direction;
                            self.rack_blocked_coordinates_text =
                                updated_rack.profile.blocked_coordinates.join(", ");
                        }
                    }
                    Err(err) => {
                        self.rack_view_status =
                            format!("Could not change rack fill direction: {}", err.message);
                        self.rack_fill_direction_editor = rack.profile.fill_direction;
                    }
                }
            }
        });
        ui.horizontal(|ui| {
            ui.label("Template");
            egui::ComboBox::from_id_salt("rack_template_combo")
                .selected_text(Self::rack_authoring_template_label(
                    self.rack_authoring_template_editor,
                ))
                .show_ui(ui, |ui| {
                    for template in [
                        RackAuthoringTemplate::BenchRows,
                        RackAuthoringTemplate::PlateColumns,
                        RackAuthoringTemplate::PlateEdgeAvoidance,
                    ] {
                        ui.selectable_value(
                            &mut self.rack_authoring_template_editor,
                            template,
                            Self::rack_authoring_template_label(template),
                        );
                    }
                });
            if ui
                .button("Apply Template")
                .on_hover_text(
                    "Apply one shared rack-authoring shortcut. Bench rows clears blocked slots and uses row-major fill; plate columns clears blocked slots and uses column-major fill; plate edge avoidance blocks the outer ring and uses column-major fill.",
                )
                .clicked()
            {
                let template = self.rack_authoring_template_editor;
                let result = self
                    .engine
                    .write()
                    .unwrap()
                    .apply(Operation::ApplyRackTemplate {
                        rack_id: rack.rack_id.clone(),
                        template,
                    });
                match result {
                    Ok(op_result) => {
                        self.rack_view_status = op_result
                            .messages
                            .first()
                            .cloned()
                            .unwrap_or_else(|| "Applied rack template".to_string());
                        self.lineage_cache_valid = false;
                        self.refresh_lineage_cache_if_needed();
                        if let Some(updated_rack) = self.current_rack_snapshot() {
                            self.rack_profile_editor_kind = updated_rack.profile.kind;
                            self.rack_authoring_template_editor =
                                Self::infer_rack_authoring_template(&updated_rack.profile);
                            self.rack_fill_direction_editor = updated_rack.profile.fill_direction;
                            self.rack_custom_profile_rows = updated_rack.profile.rows.to_string();
                            self.rack_custom_profile_columns =
                                updated_rack.profile.columns.to_string();
                            self.rack_blocked_coordinates_text =
                                updated_rack.profile.blocked_coordinates.join(", ");
                        }
                    }
                    Err(err) => {
                        self.rack_view_status =
                            format!("Could not apply rack template: {}", err.message);
                    }
                }
            }
            ui.small(
                "Templates are quick shared presets layered on top of the same rack snapshot.",
            );
        });
        if self.rack_profile_editor_kind == RackProfileKind::Custom {
            ui.horizontal(|ui| {
                ui.label("Custom rows");
                ui.add(
                    egui::TextEdit::singleline(&mut self.rack_custom_profile_rows)
                        .desired_width(56.0),
                )
                .on_hover_text(
                    "Number of rack rows for the custom A1-style profile",
                );
                ui.label("Columns");
                ui.add(
                    egui::TextEdit::singleline(&mut self.rack_custom_profile_columns)
                        .desired_width(56.0),
                )
                .on_hover_text("Number of rack columns for the custom profile");
                if ui
                    .button("Apply Custom")
                    .on_hover_text(
                        "Apply the custom rack geometry and reflow existing occupied positions under the same order-preserving rules",
                    )
                    .clicked()
                {
                    let parsed_rows = self
                        .rack_custom_profile_rows
                        .trim()
                        .parse::<usize>()
                        .map_err(|e| e.to_string());
                    let parsed_columns = self
                        .rack_custom_profile_columns
                        .trim()
                        .parse::<usize>()
                        .map_err(|e| e.to_string());
                    match (parsed_rows, parsed_columns) {
                        (Ok(rows), Ok(columns)) => {
                            let result = self
                                .engine
                                .write()
                                .unwrap()
                                .apply(Operation::SetRackProfileCustom {
                                    rack_id: rack.rack_id.clone(),
                                    rows,
                                    columns,
                                });
                            match result {
                                Ok(op_result) => {
                                    self.rack_view_status = op_result
                                        .messages
                                        .first()
                                        .cloned()
                                        .unwrap_or_else(|| {
                                            format!(
                                                "Updated rack '{}' to custom {}x{}",
                                                rack.rack_id, rows, columns
                                            )
                                        });
                                    self.lineage_cache_valid = false;
                                    self.refresh_lineage_cache_if_needed();
                                    self.rack_profile_editor_kind = RackProfileKind::Custom;
                                    if let Some(updated_rack) = self.current_rack_snapshot() {
                                        self.rack_authoring_template_editor =
                                            Self::infer_rack_authoring_template(
                                                &updated_rack.profile,
                                            );
                                        self.rack_custom_profile_rows =
                                            updated_rack.profile.rows.to_string();
                                        self.rack_custom_profile_columns =
                                            updated_rack.profile.columns.to_string();
                                        self.rack_fill_direction_editor =
                                            updated_rack.profile.fill_direction;
                                        self.rack_blocked_coordinates_text = updated_rack
                                            .profile
                                            .blocked_coordinates
                                            .join(", ");
                                    }
                                }
                                Err(err) => {
                                    self.rack_view_status = format!(
                                        "Could not apply custom rack profile: {}",
                                        err.message
                                    );
                                }
                            }
                        }
                        (Err(err), _) => {
                            self.rack_view_status =
                                format!("Could not parse custom rows: {err}");
                        }
                        (_, Err(err)) => {
                            self.rack_view_status =
                                format!("Could not parse custom columns: {err}");
                        }
                    }
                }
                ui.small("Use A1-style coordinates; row labels continue beyond Z as AA, AB, ...");
            });
        }
        ui.horizontal(|ui| {
            ui.label("Blocked");
            ui.add(
                egui::TextEdit::singleline(&mut self.rack_blocked_coordinates_text)
                    .desired_width(340.0),
            )
            .on_hover_text(
                "Comma- or whitespace-separated rack coordinates to reserve from placement, for example A1 B2 AA3",
            );
            if ui
                .button("Apply Blocked")
                .on_hover_text(
                    "Persist the blocked rack coordinates and reflow occupied positions onto the remaining available slots",
                )
                .clicked()
            {
                let blocked_coordinates =
                    Self::parse_blocked_coordinate_text(&self.rack_blocked_coordinates_text);
                let result = self
                    .engine
                    .write()
                    .unwrap()
                    .apply(Operation::SetRackBlockedCoordinates {
                        rack_id: rack.rack_id.clone(),
                        blocked_coordinates,
                    });
                match result {
                    Ok(op_result) => {
                        self.rack_view_status = op_result
                            .messages
                            .first()
                            .cloned()
                            .unwrap_or_else(|| "Updated blocked rack coordinates".to_string());
                        self.lineage_cache_valid = false;
                        self.refresh_lineage_cache_if_needed();
                        if let Some(updated_rack) = self.current_rack_snapshot() {
                            self.rack_authoring_template_editor =
                                Self::infer_rack_authoring_template(&updated_rack.profile);
                            self.rack_blocked_coordinates_text =
                                updated_rack.profile.blocked_coordinates.join(", ");
                        }
                    }
                    Err(err) => {
                        self.rack_view_status =
                            format!("Could not update blocked rack coordinates: {}", err.message);
                    }
                }
            }
            if ui
                .button("Clear Blocked")
                .on_hover_text("Clear all blocked rack coordinates for this rack")
                .clicked()
            {
                let result = self
                    .engine
                    .write()
                    .unwrap()
                    .apply(Operation::SetRackBlockedCoordinates {
                        rack_id: rack.rack_id.clone(),
                        blocked_coordinates: vec![],
                    });
                match result {
                    Ok(op_result) => {
                        self.rack_view_status = op_result
                            .messages
                            .first()
                            .cloned()
                            .unwrap_or_else(|| "Cleared blocked rack coordinates".to_string());
                        self.lineage_cache_valid = false;
                        self.refresh_lineage_cache_if_needed();
                        if let Some(updated_rack) = self.current_rack_snapshot() {
                            self.rack_authoring_template_editor =
                                Self::infer_rack_authoring_template(&updated_rack.profile);
                        }
                        self.rack_blocked_coordinates_text.clear();
                    }
                    Err(err) => {
                        self.rack_view_status =
                            format!("Could not clear blocked rack coordinates: {}", err.message);
                    }
                }
            }
        });
        let sorted_entries = Self::rack_sorted_entries(&rack);
        let blocked_coordinates = rack
            .profile
            .blocked_coordinates
            .iter()
            .cloned()
            .collect::<HashSet<_>>();
        let pointer_pos = ui.ctx().input(|i| i.pointer.hover_pos());
        let pointer_released = ui.ctx().input(|i| i.pointer.any_released());
        let escape_pressed = ui.ctx().input(|i| i.key_pressed(Key::Escape));
        let now_seconds = ui.ctx().input(|i| i.time);
        let mut pending_drag_start: Option<RackDragState> = None;
        let mut drop_target_coordinate: Option<String> = None;
        let drag_highlight_color = self.rack_view_drag_state.as_ref().map(|drag| match drag {
            RackDragState::Sample { arrangement_id, .. }
            | RackDragState::Samples { arrangement_id, .. }
            | RackDragState::ArrangementBlock { arrangement_id, .. } => {
                Self::rack_color_for_arrangement(arrangement_id)
            }
            RackDragState::ArrangementBlocks {
                arrangement_ids, ..
            } => arrangement_ids
                .first()
                .map(|arrangement_id| Self::rack_color_for_arrangement(arrangement_id))
                .unwrap_or(ui.visuals().strong_text_color()),
        });
        let previous_hover_target = self.rack_view_hover_target_coordinate.clone();
        let mut entry_by_coordinate: HashMap<String, crate::engine::RackPlacementEntry> =
            HashMap::new();
        for (_, coordinate, entry) in &sorted_entries {
            entry_by_coordinate.insert(coordinate.clone(), entry.clone());
        }
        let ghost_preview = self.rack_view_drag_state.as_ref().and_then(|drag| {
            previous_hover_target.as_deref().and_then(|target| {
                Self::rack_ghost_preview_map(&rack, &sorted_entries, drag, target)
            })
        });
        let selected_arrangement_ids = Self::rack_selected_arrangement_ids_in_order(
            &sorted_entries,
            &self.rack_view_selected_arrangement_ids,
        );
        let selected_sample_coordinates = Self::rack_selected_coordinates_in_order(
            &sorted_entries,
            &self.rack_view_selected_coordinates,
        );
        let selected_sample_arrangement_id = Self::rack_selected_samples_arrangement_id(
            &sorted_entries,
            &self.rack_view_selected_coordinates,
        );
        let preview_fill_alpha = Self::rack_drag_preview_animation_alpha(now_seconds);
        let arrangement_ids = sorted_entries
            .iter()
            .map(|(_, _, entry)| entry.arrangement_id.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        self.render_rack_help_strip(ui, &arrangement_ids);
        ui.horizontal_wrapped(|ui| {
            ui.label("Arrangement blocks");
            for arrangement_id in &arrangement_ids {
                let selected = self
                    .rack_view_selected_arrangement_ids
                    .contains(arrangement_id);
                let color = Self::rack_color_for_arrangement(arrangement_id);
                let button = egui::Button::new(
                    egui::RichText::new(arrangement_id.clone()).color(color).strong(),
                )
                .selected(selected)
                .sense(egui::Sense::click_and_drag());
                let response = ui
                    .add(button)
                    .on_hover_text("Drag this arrangement block onto another rack slot, or click once to select it for a later target click. Command/Ctrl-click toggles multi-selection.");
                if response.drag_started() {
                    if selected && selected_arrangement_ids.len() > 1 {
                        if let Some(first_selected) = selected_arrangement_ids.first()
                            && let Some(from_coordinate) = Self::rack_first_coordinate_for_arrangement(
                                &sorted_entries,
                                first_selected,
                            ) {
                                pending_drag_start = Some(RackDragState::ArrangementBlocks {
                                    arrangement_ids: selected_arrangement_ids.clone(),
                                    from_coordinate,
                                });
                            }
                    } else if let Some(from_coordinate) =
                        Self::rack_first_coordinate_for_arrangement(&sorted_entries, arrangement_id)
                    {
                        pending_drag_start = Some(RackDragState::ArrangementBlock {
                            arrangement_id: arrangement_id.clone(),
                            from_coordinate,
                        });
                    }
                } else if response.clicked() {
                    let additive = ui.ctx().input(|i| i.modifiers.command || i.modifiers.ctrl);
                    if additive {
                        if selected {
                            self.rack_view_selected_arrangement_ids.remove(arrangement_id);
                        } else {
                            self.rack_view_selected_arrangement_ids
                                .insert(arrangement_id.clone());
                        }
                    } else {
                        if selected && self.rack_view_selected_arrangement_ids.len() == 1 {
                            self.rack_view_selected_arrangement_ids.clear();
                        } else {
                            self.rack_view_selected_arrangement_ids.clear();
                            self.rack_view_selected_arrangement_ids
                                .insert(arrangement_id.clone());
                        }
                    }
                    self.rack_view_selected_coordinates.clear();
                }
            }
            if ui
                .button("Clear selection")
                .on_hover_text("Clear current sample/block selection or cancel an in-progress rack drag")
                .clicked()
            {
                self.rack_view_selected_coordinates.clear();
                self.rack_view_selected_arrangement_ids.clear();
                self.rack_view_drag_state = None;
            }
        });

        if !self.rack_view_status.trim().is_empty() {
            ui.small(self.rack_view_status.clone());
        }
        ui.separator();
        let rack_scroll_output = egui::ScrollArea::both()
            .id_salt("rack_grid_scroll")
            .auto_shrink([false, false])
            .scroll_offset(self.rack_view_scroll_offset)
            .show(ui, |ui| {
            egui::Grid::new("rack_grid")
                .spacing(egui::vec2(6.0, 6.0))
                .show(ui, |ui| {
                    ui.label("");
                    for column in 0..rack.profile.columns {
                        ui.monospace((column + 1).to_string());
                    }
                    ui.end_row();
                    for row in 0..rack.profile.rows {
                        ui.monospace(GentleEngine::rack_row_label_from_index(row));
                        for column in 0..rack.profile.columns {
                            let coordinate =
                                Self::rack_coordinate_for_slot(&rack.profile, row, column);
                            if blocked_coordinates.contains(&coordinate) {
                                ui.add_sized(
                                    egui::vec2(96.0, 44.0),
                                    egui::Button::new(
                                        egui::RichText::new(format!("Blocked\n{coordinate}"))
                                            .weak()
                                            .small(),
                                    )
                                    .sense(egui::Sense::hover()),
                                )
                                .on_hover_text(
                                    "This rack position is blocked/reserved and is excluded from placement order",
                                );
                                continue;
                            }

                            let preview_cell =
                                ghost_preview.as_ref().and_then(|map| map.get(&coordinate));
                            let current_entry = entry_by_coordinate.get(&coordinate).cloned();
                            let effective_entry = preview_cell
                                .and_then(|cell| cell.predicted_entry.clone())
                                .or(current_entry.clone());

                            if let Some(entry) = effective_entry {
                                let arrangement_color =
                                    Self::rack_color_for_arrangement(&entry.arrangement_id);
                                let display = Self::rack_entry_display(&entry);
                                let selected =
                                    self.rack_view_selected_coordinates.contains(&coordinate);
                                let mut button = egui::Button::new(
                                    egui::RichText::new(display)
                                        .color(arrangement_color)
                                        .small(),
                                )
                                .selected(selected)
                                .sense(egui::Sense::click_and_drag());
                                if preview_cell.is_some() {
                                    button =
                                        button.fill(arrangement_color.gamma_multiply(preview_fill_alpha));
                                }
                                let hover_text = if preview_cell.is_some() {
                                    "Ghost preview after drop. Drag this sample to another rack slot, or click once to select it for a later target click"
                                } else {
                                    "Drag this sample to another rack slot, or click once to select it for a later target click"
                                };
                                let response = ui
                                    .add_sized(egui::vec2(96.0, 44.0), button)
                                    .on_hover_text(hover_text);
                                let is_drop_target = self.rack_view_drag_state.is_some()
                                    && pointer_pos
                                        .map(|pos| response.rect.contains(pos))
                                        .unwrap_or(false);
                                if is_drop_target {
                                    drop_target_coordinate = Some(coordinate.clone());
                                    ui.painter().rect_stroke(
                                        response.rect.expand(1.5),
                                        6.0,
                                        egui::Stroke::new(
                                            2.0,
                                            drag_highlight_color.unwrap_or(arrangement_color),
                                        ),
                                        egui::StrokeKind::Outside,
                                    );
                                } else if preview_cell.is_some() {
                                    ui.painter().rect_stroke(
                                        response.rect.expand(1.0),
                                        6.0,
                                        egui::Stroke::new(
                                            1.5,
                                            arrangement_color.gamma_multiply(0.85),
                                        ),
                                        egui::StrokeKind::Outside,
                                    );
                                }
                                if let Some(drop_ghost) = self.rack_view_recent_drop_ghost.as_ref()
                                    && drop_ghost.rack_id == rack.rack_id
                                        && drop_ghost.changed_coordinates.contains(&coordinate)
                                    {
                                        let fade = 1.0
                                            - (drop_ghost.started_at.elapsed().as_secs_f32() / 0.7)
                                                .clamp(0.0, 1.0);
                                        if fade > 0.0 {
                                            ui.painter().rect_filled(
                                                response.rect.shrink(3.0),
                                                5.0,
                                                arrangement_color.gamma_multiply(0.05 + 0.10 * fade),
                                            );
                                            ui.painter().rect_stroke(
                                                response.rect.expand(1.0),
                                                6.0,
                                                egui::Stroke::new(
                                                    1.5,
                                                    arrangement_color.gamma_multiply(0.35 + 0.45 * fade),
                                                ),
                                                egui::StrokeKind::Outside,
                                            );
                                            ui.ctx().request_repaint();
                                        }
                                    }
                                if current_entry.is_some() && response.drag_started() {
                                    let current_arrangement_id = current_entry
                                        .as_ref()
                                        .map(|value| value.arrangement_id.clone())
                                        .unwrap_or_else(|| entry.arrangement_id.clone());
                                    if selected
                                        && selected_sample_coordinates.len() > 1
                                        && selected_sample_arrangement_id
                                            .as_deref()
                                            .map(|value| value == current_arrangement_id)
                                            .unwrap_or(false)
                                    {
                                        pending_drag_start = Some(RackDragState::Samples {
                                            arrangement_id: current_arrangement_id,
                                            from_coordinates: selected_sample_coordinates.clone(),
                                        });
                                    } else {
                                        pending_drag_start = Some(RackDragState::Sample {
                                            from_coordinate: coordinate.clone(),
                                            arrangement_id: current_entry
                                                .as_ref()
                                                .map(|value| value.arrangement_id.clone())
                                                .unwrap_or_else(|| entry.arrangement_id.clone()),
                                            role_label: current_entry
                                                .as_ref()
                                                .map(|value| value.role_label.clone())
                                                .unwrap_or_else(|| entry.role_label.clone()),
                                        });
                                    }
                                } else if response.clicked() && self.rack_view_drag_state.is_none() {
                                    let additive = ui.ctx().input(|i| i.modifiers.command || i.modifiers.ctrl);
                                    if additive {
                                        self.rack_view_selected_arrangement_ids.clear();
                                        let current_arrangement_id = current_entry
                                            .as_ref()
                                            .map(|value| value.arrangement_id.clone())
                                            .unwrap_or_else(|| entry.arrangement_id.clone());
                                        if selected {
                                            self.rack_view_selected_coordinates.remove(&coordinate);
                                        } else if let Some(existing_arrangement_id) =
                                            selected_sample_arrangement_id.as_deref()
                                        {
                                            if existing_arrangement_id != current_arrangement_id {
                                                self.rack_view_selected_coordinates.clear();
                                            }
                                            self.rack_view_selected_coordinates
                                                .insert(coordinate.clone());
                                        } else {
                                            self.rack_view_selected_coordinates
                                                .insert(coordinate.clone());
                                        }
                                    } else if !selected_arrangement_ids.is_empty() {
                                        if selected_arrangement_ids.len() == 1 {
                                            if let Some(from_coordinate) = Self::rack_first_coordinate_for_arrangement(
                                                &sorted_entries,
                                                &selected_arrangement_ids[0],
                                            ) {
                                                self.apply_rack_move(
                                                    &rack.rack_id,
                                                    &from_coordinate,
                                                    &coordinate,
                                                    true,
                                                );
                                            }
                                        } else {
                                            self.apply_rack_move_blocks(
                                                &rack.rack_id,
                                                &selected_arrangement_ids,
                                                &coordinate,
                                            );
                                        }
                                        self.rack_view_selected_arrangement_ids.clear();
                                    } else if !selected_sample_coordinates.is_empty() {
                                        if selected_sample_coordinates.len() == 1
                                            && selected_sample_coordinates[0] == coordinate
                                        {
                                            self.rack_view_selected_coordinates.clear();
                                        } else if selected_sample_coordinates.len() == 1 {
                                            self.apply_rack_move(
                                                &rack.rack_id,
                                                &selected_sample_coordinates[0],
                                                &coordinate,
                                                false,
                                            );
                                            self.rack_view_selected_coordinates.clear();
                                        } else {
                                            self.apply_rack_move_samples(
                                                &rack.rack_id,
                                                &selected_sample_coordinates,
                                                &coordinate,
                                            );
                                            self.rack_view_selected_coordinates.clear();
                                        }
                                    } else if current_entry.is_some() {
                                        self.rack_view_selected_coordinates.clear();
                                        self.rack_view_selected_coordinates
                                            .insert(coordinate.clone());
                                        self.rack_view_selected_arrangement_ids.clear();
                                    }
                                }
                            } else {
                                let mut button =
                                    egui::Button::new(egui::RichText::new(coordinate.clone()).weak());
                                if preview_cell.is_some() {
                                    button = button.fill(ui.visuals().faint_bg_color);
                                }
                                let response = ui
                                    .add_sized(egui::vec2(96.0, 44.0), button)
                                    .on_hover_text(if preview_cell.is_some() {
                                        "Ghost preview after drop. Drop a dragged sample or arrangement block here, or click as the target for an already selected move"
                                    } else {
                                        "Drop a dragged sample or arrangement block here, or click as the target for an already selected move"
                                    });
                                let is_drop_target = self.rack_view_drag_state.is_some()
                                    && pointer_pos
                                        .map(|pos| response.rect.contains(pos))
                                        .unwrap_or(false);
                                if is_drop_target {
                                    drop_target_coordinate = Some(coordinate.clone());
                                    ui.painter().rect_stroke(
                                        response.rect.expand(1.5),
                                        6.0,
                                        egui::Stroke::new(
                                            2.0,
                                            drag_highlight_color
                                                .unwrap_or(ui.visuals().strong_text_color()),
                                        ),
                                        egui::StrokeKind::Outside,
                                    );
                                } else if preview_cell.is_some() {
                                    ui.painter().rect_stroke(
                                        response.rect.expand(1.0),
                                        6.0,
                                        egui::Stroke::new(
                                            1.5,
                                            drag_highlight_color
                                                .unwrap_or(ui.visuals().strong_text_color())
                                                .gamma_multiply(0.8),
                                        ),
                                        egui::StrokeKind::Outside,
                                    );
                                }
                                if let Some(drop_ghost) = self.rack_view_recent_drop_ghost.as_ref()
                                    && drop_ghost.rack_id == rack.rack_id
                                        && drop_ghost.changed_coordinates.contains(&coordinate)
                                    {
                                        let fade = 1.0
                                            - (drop_ghost.started_at.elapsed().as_secs_f32() / 0.7)
                                                .clamp(0.0, 1.0);
                                        if fade > 0.0 {
                                            ui.painter().rect_stroke(
                                                response.rect.expand(1.0),
                                                6.0,
                                                egui::Stroke::new(
                                                    1.2,
                                                    drag_highlight_color
                                                        .unwrap_or(ui.visuals().strong_text_color())
                                                        .gamma_multiply(0.30 + 0.40 * fade),
                                                ),
                                                egui::StrokeKind::Outside,
                                            );
                                            ui.ctx().request_repaint();
                                        }
                                    }
                                if response.clicked() && self.rack_view_drag_state.is_none() {
                                    if !selected_arrangement_ids.is_empty() {
                                        if selected_arrangement_ids.len() == 1 {
                                            if let Some(from_coordinate) = Self::rack_first_coordinate_for_arrangement(
                                                &sorted_entries,
                                                &selected_arrangement_ids[0],
                                            ) {
                                                self.apply_rack_move(
                                                    &rack.rack_id,
                                                    &from_coordinate,
                                                    &coordinate,
                                                    true,
                                                );
                                            }
                                        } else {
                                            self.apply_rack_move_blocks(
                                                &rack.rack_id,
                                                &selected_arrangement_ids,
                                                &coordinate,
                                            );
                                        }
                                        self.rack_view_selected_arrangement_ids.clear();
                                    } else if !selected_sample_coordinates.is_empty() {
                                        if selected_sample_coordinates.len() == 1 {
                                            self.apply_rack_move(
                                                &rack.rack_id,
                                                &selected_sample_coordinates[0],
                                                &coordinate,
                                                false,
                                            );
                                        } else {
                                            self.apply_rack_move_samples(
                                                &rack.rack_id,
                                                &selected_sample_coordinates,
                                                &coordinate,
                                            );
                                        }
                                        self.rack_view_selected_coordinates.clear();
                                    }
                                }
                            }
                        }
                        ui.end_row();
                    }
                });
        });
        let max_scroll_x =
            (rack_scroll_output.content_size.x - rack_scroll_output.inner_rect.width()).max(0.0);
        let max_scroll_y =
            (rack_scroll_output.content_size.y - rack_scroll_output.inner_rect.height()).max(0.0);
        let mut next_rack_scroll_offset = rack_scroll_output.state.offset;
        next_rack_scroll_offset.x = next_rack_scroll_offset.x.clamp(0.0, max_scroll_x);
        next_rack_scroll_offset.y = next_rack_scroll_offset.y.clamp(0.0, max_scroll_y);
        if self.rack_view_drag_state.is_some()
            && let Some(pointer) = pointer_pos
        {
            let autoscrolled = Self::rack_autoscroll_offset_for_pointer(
                next_rack_scroll_offset,
                pointer,
                rack_scroll_output.inner_rect,
                max_scroll_x,
                max_scroll_y,
            );
            if (autoscrolled - next_rack_scroll_offset).length_sq() > 0.01 {
                next_rack_scroll_offset = autoscrolled;
                ui.ctx().request_repaint();
            }
        }
        self.rack_view_scroll_offset = next_rack_scroll_offset;
        if let Some(drag) = pending_drag_start {
            self.rack_view_drag_state = Some(drag);
            self.rack_view_selected_coordinates.clear();
        }
        if escape_pressed && self.rack_view_drag_state.is_some() {
            self.rack_view_drag_state = None;
            self.rack_view_status = "Canceled rack drag.".to_string();
        }
        if let Some(drag) = self.rack_view_drag_state.as_ref() {
            ui.small(Self::rack_drag_status_text(
                drag,
                drop_target_coordinate.as_deref(),
            ));
        }
        if pointer_released && let Some(drag) = self.rack_view_drag_state.clone() {
            match drag {
                RackDragState::Sample {
                    from_coordinate, ..
                } => {
                    if let Some(target) = drop_target_coordinate.as_deref() {
                        self.apply_rack_move(&rack.rack_id, &from_coordinate, target, false);
                    }
                }
                RackDragState::Samples {
                    from_coordinates, ..
                } => {
                    if let Some(target) = drop_target_coordinate.as_deref() {
                        self.apply_rack_move_samples(&rack.rack_id, &from_coordinates, target);
                    }
                }
                RackDragState::ArrangementBlock {
                    from_coordinate, ..
                } => {
                    if let Some(target) = drop_target_coordinate.as_deref() {
                        self.apply_rack_move(&rack.rack_id, &from_coordinate, target, true);
                    }
                }
                RackDragState::ArrangementBlocks {
                    arrangement_ids, ..
                } => {
                    if let Some(target) = drop_target_coordinate.as_deref() {
                        self.apply_rack_move_blocks(&rack.rack_id, &arrangement_ids, target);
                    }
                }
            }
            self.rack_view_drag_state = None;
        }
        if self.rack_view_drag_state.is_some() {
            if previous_hover_target != drop_target_coordinate {
                ui.ctx().request_repaint();
            }
            self.rack_view_hover_target_coordinate = drop_target_coordinate.clone();
        } else {
            self.rack_view_hover_target_coordinate = None;
        }
        if let Some(drop_ghost) = self.rack_view_recent_drop_ghost.as_ref()
            && (drop_ghost.rack_id != rack.rack_id
                || drop_ghost.started_at.elapsed().as_secs_f32() >= 0.7)
        {
            self.rack_view_recent_drop_ghost = None;
        }
        ui.separator();
        ui.horizontal(|ui| {
            ui.label("Label preset");
            egui::ComboBox::from_id_salt("rack_label_sheet_preset_combo")
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
            if ui
                .button("Preview Labels...")
                .on_hover_text("Open an in-app preview of the deterministic label sheet for the currently shown rack before exporting it")
                .clicked()
            {
                self.open_rack_labels_preview_dialog(&rack.rack_id, None);
            }
            ui.separator();
            ui.label("Physical template");
            egui::ComboBox::from_id_salt("rack_physical_template_combo")
                .selected_text(Self::rack_physical_template_label(
                    self.rack_physical_template_kind,
                ))
                .show_ui(ui, |ui| {
                    for template in [
                        RackPhysicalTemplateKind::StoragePcrTubeRack,
                        RackPhysicalTemplateKind::PipettingPcrTubeRack,
                        RackPhysicalTemplateKind::CellCulturePlate,
                    ] {
                        ui.selectable_value(
                            &mut self.rack_physical_template_kind,
                            template,
                            Self::rack_physical_template_label(template),
                        );
                    }
                });
            ui.label("Carrier preset");
            egui::ComboBox::from_id_salt("rack_carrier_label_preset_combo")
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
            if ui
                .button("Fabrication SVG...")
                .on_hover_text("Export a top-view fabrication/planning SVG for the current rack using the selected physical carrier template")
                .clicked()
            {
                self.prompt_export_rack_fabrication_svg(&rack.rack_id);
            }
            if ui
                .button("Isometric SVG...")
                .on_hover_text("Export a presentation-grade pseudo-3D rack SVG for the current rack using the selected physical carrier template")
                .clicked()
            {
                self.prompt_export_rack_isometric_svg(&rack.rack_id);
            }
            if ui
                .button("OpenSCAD...")
                .on_hover_text("Export one parameterized OpenSCAD file for the current rack using the selected physical carrier template")
                .clicked()
            {
                self.prompt_export_rack_openscad(&rack.rack_id);
            }
            if ui
                .button("Carrier labels SVG...")
                .on_hover_text("Export one carrier-matched front-strip plus module-label SVG sheet for the current rack using the selected physical carrier template")
                .clicked()
            {
                self.prompt_export_rack_carrier_labels_svg(&rack.rack_id);
            }
            if ui
                .button("Simulation JSON...")
                .on_hover_text("Export one machine-readable physical rack geometry/placement JSON for downstream simulation adapters using the selected physical carrier template")
                .clicked()
            {
                self.prompt_export_rack_simulation_json(&rack.rack_id);
            }
        });
        close_requested
    }

    pub(super) fn render_rack_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_rack_dialog {
            return;
        }
        let mut open = self.show_rack_dialog;
        let title = if self.rack_view_rack_id.trim().is_empty() {
            "Rack".to_string()
        } else {
            format!("Rack — {}", self.rack_view_rack_id.trim())
        };
        let viewport_id = Self::rack_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            title.clone(),
            egui::Id::new(("hosted_rack_window", viewport_id)),
            viewport_id,
            Vec2::new(1000.0, 760.0),
            Vec2::new(760.0, 520.0),
        );
        if ctx.embed_viewports() {
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                egui::ScrollArea::both()
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        close_requested = self.render_rack_contents(ui);
                    });
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            if close_requested || ctx.input(|i| i.key_pressed(Key::Escape)) {
                open = false;
            }
            self.show_rack_dialog = open;
            return;
        }
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                    egui::ScrollArea::both()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            close_requested = self.render_rack_contents(ui);
                        });
                });
                if close_requested {
                    open = false;
                }
            } else {
                let mut close_requested = false;
                crate::egui_compat::show_central_panel(
                    &mut *ctx,
                    egui::CentralPanel::default(),
                    |ui| {
                        close_requested = self.render_rack_contents(ui);
                    },
                );
                if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            open = false;
        }
        self.show_rack_dialog = open;
    }

    pub(super) fn open_place_arrangement_on_rack_dialog(&mut self, arrangement_id: &str) {
        self.place_arrangement_source_id = arrangement_id.trim().to_string();
        if self.place_arrangement_target_rack_id.trim().is_empty()
            && let Some(first_rack) = self.lineage_racks.first()
        {
            self.place_arrangement_target_rack_id = first_rack.rack_id.clone();
        }
        self.place_arrangement_status.clear();
        self.show_place_arrangement_rack_dialog = true;
    }

    pub(super) fn submit_place_arrangement_on_rack(&mut self) {
        let arrangement_id = self.place_arrangement_source_id.trim().to_string();
        let rack_id = self.place_arrangement_target_rack_id.trim().to_string();
        if arrangement_id.is_empty() || rack_id.is_empty() {
            self.place_arrangement_status =
                "Choose both an arrangement and a target rack.".to_string();
            return;
        }
        let result = {
            self.engine
                .write()
                .unwrap()
                .apply(Operation::PlaceArrangementOnRack {
                    arrangement_id: arrangement_id.clone(),
                    rack_id: rack_id.clone(),
                })
        };
        match result {
            Ok(op_result) => {
                self.place_arrangement_status =
                    op_result.messages.first().cloned().unwrap_or_else(|| {
                        format!("Placed arrangement '{}' onto '{}'", arrangement_id, rack_id)
                    });
                self.lineage_cache_valid = false;
                self.refresh_lineage_cache_if_needed();
                self.show_place_arrangement_rack_dialog = false;
            }
            Err(err) => {
                self.place_arrangement_status = format!(
                    "Could not place arrangement '{}' onto rack '{}': {}",
                    arrangement_id, rack_id, err.message
                );
            }
        }
    }

    pub(super) fn render_place_arrangement_on_rack_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_place_arrangement_rack_dialog {
            return;
        }
        let mut open = self.show_place_arrangement_rack_dialog;
        let mut close_requested = false;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Place Arrangement on Existing Rack",
            egui::Id::new("place_arrangement_on_rack_modal"),
        )
        .default_size(Vec2::new(560.0, 220.0));
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            ui.label("Append one arrangement as a contiguous block onto an existing rack. Existing blocks stay intact.");
            ui.horizontal(|ui| {
                ui.label("Arrangement");
                ui.monospace(self.place_arrangement_source_id.clone());
            });
            ui.horizontal(|ui| {
                ui.label("Target rack");
                egui::ComboBox::from_id_salt("place_arrangement_rack_combo")
                    .selected_text(if self.place_arrangement_target_rack_id.trim().is_empty() {
                        "(choose rack)".to_string()
                    } else {
                        self.place_arrangement_target_rack_id.clone()
                    })
                    .show_ui(ui, |ui| {
                        for rack in &self.lineage_racks {
                            ui.selectable_value(
                                &mut self.place_arrangement_target_rack_id,
                                rack.rack_id.clone(),
                                format!(
                                    "{} — {} ({}, {} occupied, {} arrangement block(s))",
                                    rack.rack_id,
                                    rack.name,
                                    rack.profile,
                                    rack.occupied_positions,
                                    rack.arrangement_ids.len()
                                ),
                            );
                        }
                    });
            });
            if !self.place_arrangement_status.trim().is_empty() {
                ui.small(self.place_arrangement_status.clone());
            }
            ui.separator();
            ui.horizontal(|ui| {
                if ui.button("Place").clicked() {
                    self.submit_place_arrangement_on_rack();
                }
                if ui.button("Cancel").clicked() {
                    close_requested = true;
                }
            });
        });
        if close_requested {
            open = false;
        }
        self.show_place_arrangement_rack_dialog = open;
    }
}
