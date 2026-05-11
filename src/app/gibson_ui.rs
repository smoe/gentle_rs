//! Gibson specialist GUI state and rendering helpers.
//!
//! This module is a move-only extraction from `app.rs`: it keeps the
//! Gibson specialist, preview/apply/export handlers, and related
//! UI-specific state close to `GENtleApp` while reducing the top-level
//! app monolith.

use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum GibsonUiOpeningMode {
    ExistingTermini,
    DefinedSite,
}

impl GibsonUiOpeningMode {
    fn as_plan_mode(self) -> &'static str {
        match self {
            Self::ExistingTermini => "existing_termini",
            Self::DefinedSite => "defined_site",
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::ExistingTermini => "Use existing termini",
            Self::DefinedSite => "Defined cut/opening",
        }
    }

    fn from_plan_mode(raw: &str) -> Self {
        match raw.trim().to_ascii_lowercase().as_str() {
            "existing_termini" | "existing-termini" => Self::ExistingTermini,
            _ => Self::DefinedSite,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum GibsonUiInsertOrientation {
    Forward,
    Reverse,
}

impl GibsonUiInsertOrientation {
    fn as_plan_orientation(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Forward => "Forward",
            Self::Reverse => "Reverse",
        }
    }

    fn from_plan_orientation(raw: &str) -> Self {
        match raw.trim().to_ascii_lowercase().as_str() {
            "reverse" | "rev" | "reverse_complement" => Self::Reverse,
            _ => Self::Forward,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct GibsonUiInsertRow {
    pub(super) seq_id: String,
    pub(super) orientation: GibsonUiInsertOrientation,
}

impl Default for GibsonUiInsertRow {
    fn default() -> Self {
        Self {
            seq_id: String::new(),
            orientation: GibsonUiInsertOrientation::Forward,
        }
    }
}

impl GENtleApp {
    pub(super) fn build_gibson_plan_from_ui(
        &self,
    ) -> std::result::Result<GibsonAssemblyPlan, String> {
        fn parse_usize_field(label: &str, raw: &str) -> std::result::Result<usize, String> {
            raw.trim()
                .parse::<usize>()
                .map_err(|e| format!("Could not parse {label} '{raw}': {e}"))
        }

        fn parse_f64_field(label: &str, raw: &str) -> std::result::Result<f64, String> {
            raw.trim()
                .parse::<f64>()
                .map_err(|e| format!("Could not parse {label} '{raw}': {e}"))
        }

        let destination_seq_id = self.gibson_destination_seq_id.trim();
        if destination_seq_id.is_empty() {
            return Err("Choose a destination sequence for Gibson planning".to_string());
        }
        let insert_rows = self.gibson_ui_insert_rows();
        if insert_rows.is_empty() {
            return Err("Choose at least one insert sequence for Gibson planning".to_string());
        }
        if insert_rows
            .iter()
            .any(|row| destination_seq_id == row.seq_id.trim())
        {
            return Err(
                "Destination and insert sequence ids must differ for Gibson planning".to_string(),
            );
        }

        let destination_topology = {
            let engine = self.engine.read().expect("Engine lock poisoned");
            let destination = engine
                .state()
                .sequences
                .get(destination_seq_id)
                .ok_or_else(|| {
                    format!("Destination sequence '{destination_seq_id}' does not exist in state")
                })?;
            if destination.is_circular() {
                "circular".to_string()
            } else {
                "linear".to_string()
            }
        };

        let (opening_start, opening_end) = match self.gibson_opening_mode {
            GibsonUiOpeningMode::ExistingTermini => (None, None),
            GibsonUiOpeningMode::DefinedSite => (
                Some(parse_usize_field(
                    "left cut edge (0-based)",
                    &self.gibson_opening_start_0based,
                )?),
                Some(parse_usize_field(
                    "right cut edge (0-based)",
                    &self.gibson_opening_end_0based_exclusive,
                )?),
            ),
        };

        let overlap_bp_min = parse_usize_field("overlap_bp_min", &self.gibson_overlap_bp_min)?;
        let overlap_bp_max = parse_usize_field("overlap_bp_max", &self.gibson_overlap_bp_max)?;
        let priming_length_min_bp = parse_usize_field(
            "priming_segment_min_length_bp",
            &self.gibson_priming_length_min_bp,
        )?;
        let priming_length_max_bp = parse_usize_field(
            "priming_segment_max_length_bp",
            &self.gibson_priming_length_max_bp,
        )?;
        if overlap_bp_min > overlap_bp_max {
            return Err(format!(
                "overlap bp minimum {} exceeds maximum {}",
                overlap_bp_min, overlap_bp_max
            ));
        }
        if priming_length_min_bp > priming_length_max_bp {
            return Err(format!(
                "priming length minimum {} exceeds maximum {}",
                priming_length_min_bp, priming_length_max_bp
            ));
        }
        let minimum_overlap_tm_celsius = parse_f64_field(
            "minimum_overlap_tm_celsius",
            &self.gibson_minimum_overlap_tm_celsius,
        )?;
        let priming_tm_min_celsius = parse_f64_field(
            "priming_segment_tm_min_celsius",
            &self.gibson_priming_tm_min_celsius,
        )?;
        let priming_tm_max_celsius = parse_f64_field(
            "priming_segment_tm_max_celsius",
            &self.gibson_priming_tm_max_celsius,
        )?;
        if priming_tm_min_celsius > priming_tm_max_celsius {
            return Err(format!(
                "priming Tm minimum {:.1} exceeds maximum {:.1}",
                priming_tm_min_celsius, priming_tm_max_celsius
            ));
        }
        let output_id_hint = if self.gibson_output_id_hint.trim().is_empty() {
            format!(
                "{destination_seq_id}_with_{}",
                insert_rows
                    .iter()
                    .map(|row| row.seq_id.trim())
                    .collect::<Vec<_>>()
                    .join("_")
            )
        } else {
            self.gibson_output_id_hint.trim().to_string()
        };
        let desired_unique_restriction_site_enzyme_name = (!self
            .gibson_unique_restriction_site_enzyme_name
            .trim()
            .is_empty())
        .then(|| {
            self.gibson_unique_restriction_site_enzyme_name
                .trim()
                .to_string()
        });

        let fragments = insert_rows
            .iter()
            .enumerate()
            .map(|(idx, row)| GibsonPlanFragment {
                id: format!("insert_{}", idx + 1),
                seq_id: row.seq_id.trim().to_string(),
                role: "insert".to_string(),
                orientation: row.orientation.as_plan_orientation().to_string(),
                left_end_strategy: Some(GibsonPlanEndStrategy {
                    mode: "primer_added_overlap".to_string(),
                    target_junction_id: if idx == 0 {
                        "junction_left".to_string()
                    } else {
                        format!("junction_{}_{}", idx, idx + 1)
                    },
                }),
                right_end_strategy: Some(GibsonPlanEndStrategy {
                    mode: "primer_added_overlap".to_string(),
                    target_junction_id: if idx + 1 == insert_rows.len() {
                        "junction_right".to_string()
                    } else {
                        format!("junction_{}_{}", idx + 1, idx + 2)
                    },
                }),
                source_span_1based: None,
            })
            .collect::<Vec<_>>();
        let mut assembly_order = vec![GibsonPlanAssemblyMember {
            kind: "destination_end".to_string(),
            id: "dest_left".to_string(),
        }];
        assembly_order.extend((0..insert_rows.len()).map(|idx| GibsonPlanAssemblyMember {
            kind: "fragment".to_string(),
            id: format!("insert_{}", idx + 1),
        }));
        assembly_order.push(GibsonPlanAssemblyMember {
            kind: "destination_end".to_string(),
            id: "dest_right".to_string(),
        });
        let mut junctions = vec![GibsonPlanJunction {
            id: "junction_left".to_string(),
            left_member: GibsonPlanAssemblyMember {
                kind: "destination_end".to_string(),
                id: "dest_left".to_string(),
            },
            right_member: GibsonPlanAssemblyMember {
                kind: "fragment".to_string(),
                id: "insert_1".to_string(),
            },
            required_overlap_bp: None,
            overlap_partition: Some(GibsonPlanOverlapPartition {
                left_member_bp: 0,
                right_member_bp: 0,
            }),
            overlap_source: "derive_from_destination_left_flank".to_string(),
            distinct_from: vec!["junction_right".to_string()],
        }];
        for idx in 0..insert_rows.len().saturating_sub(1) {
            junctions.push(GibsonPlanJunction {
                id: format!("junction_{}_{}", idx + 1, idx + 2),
                left_member: GibsonPlanAssemblyMember {
                    kind: "fragment".to_string(),
                    id: format!("insert_{}", idx + 1),
                },
                right_member: GibsonPlanAssemblyMember {
                    kind: "fragment".to_string(),
                    id: format!("insert_{}", idx + 2),
                },
                required_overlap_bp: None,
                overlap_partition: Some(GibsonPlanOverlapPartition {
                    left_member_bp: 0,
                    right_member_bp: 0,
                }),
                overlap_source: "designed_bridge_sequence".to_string(),
                distinct_from: vec![],
            });
        }
        junctions.push(GibsonPlanJunction {
            id: "junction_right".to_string(),
            left_member: GibsonPlanAssemblyMember {
                kind: "fragment".to_string(),
                id: format!("insert_{}", insert_rows.len()),
            },
            right_member: GibsonPlanAssemblyMember {
                kind: "destination_end".to_string(),
                id: "dest_right".to_string(),
            },
            required_overlap_bp: None,
            overlap_partition: Some(GibsonPlanOverlapPartition {
                left_member_bp: 0,
                right_member_bp: 0,
            }),
            overlap_source: "derive_from_destination_right_flank".to_string(),
            distinct_from: vec!["junction_left".to_string()],
        });

        Ok(GibsonAssemblyPlan {
            schema: "gentle.gibson_assembly_plan.v1".to_string(),
            id: format!(
                "gibson_plan_{}_{}",
                destination_seq_id,
                insert_rows
                    .iter()
                    .map(|row| row.seq_id.trim())
                    .collect::<Vec<_>>()
                    .join("_")
            ),
            title: format!(
                "Gibson assembly plan: {} + {}",
                destination_seq_id,
                insert_rows
                    .iter()
                    .map(|row| row.seq_id.trim())
                    .collect::<Vec<_>>()
                    .join(" + ")
            ),
            summary: if insert_rows.len() == 1 {
                "Destination-first single-insert Gibson preview from specialist UI".to_string()
            } else {
                format!(
                    "Destination-first multi-insert Gibson preview from specialist UI ({} inserts)",
                    insert_rows.len()
                )
            },
            destination: GibsonPlanDestination {
                seq_id: destination_seq_id.to_string(),
                topology_before_opening: destination_topology.clone(),
                opening: GibsonPlanOpening {
                    mode: self.gibson_opening_mode.as_plan_mode().to_string(),
                    label: match self.gibson_opening_mode {
                        GibsonUiOpeningMode::ExistingTermini => "existing termini".to_string(),
                        GibsonUiOpeningMode::DefinedSite => "selected opening".to_string(),
                    },
                    start_0based: opening_start,
                    end_0based_exclusive: opening_end,
                    left_end_id: "dest_left".to_string(),
                    right_end_id: "dest_right".to_string(),
                    uniqueness_requirement: "must_be_unambiguous".to_string(),
                },
            },
            product: GibsonPlanProduct {
                topology: destination_topology,
                output_id_hint,
            },
            fragments,
            assembly_order,
            junctions,
            validation_policy: GibsonPlanValidationPolicy {
                require_unambiguous_destination_opening: true,
                require_distinct_terminal_junctions: true,
                adjacency_overlap_mismatch: "error".to_string(),
                design_targets: GibsonPlanDesignTargets {
                    overlap_bp_min,
                    overlap_bp_max,
                    minimum_overlap_tm_celsius,
                    priming_segment_tm_min_celsius: priming_tm_min_celsius,
                    priming_segment_tm_max_celsius: priming_tm_max_celsius,
                    priming_segment_min_length_bp: priming_length_min_bp,
                    priming_segment_max_length_bp: priming_length_max_bp,
                    max_anneal_hits: 4,
                },
                uniqueness_checks: GibsonPlanUniquenessChecks {
                    destination_context: "warn".to_string(),
                    participating_fragments: "warn".to_string(),
                    reference_contexts: vec![],
                },
                desired_unique_restriction_site_enzyme_name,
            },
            derived_design: None,
        })
    }

    pub(super) fn load_gibson_plan_into_ui(&mut self, plan: &GibsonAssemblyPlan) {
        self.gibson_destination_seq_id = plan.destination.seq_id.trim().to_string();
        self.gibson_opening_mode =
            GibsonUiOpeningMode::from_plan_mode(&plan.destination.opening.mode);
        self.gibson_opening_start_0based = plan
            .destination
            .opening
            .start_0based
            .map(|value| value.to_string())
            .unwrap_or_default();
        self.gibson_opening_end_0based_exclusive = plan
            .destination
            .opening
            .end_0based_exclusive
            .map(|value| value.to_string())
            .unwrap_or_default();
        self.gibson_extra_inserts.clear();
        if let Some(fragment) = plan.fragments.first() {
            self.gibson_insert_seq_id = fragment.seq_id.trim().to_string();
            self.gibson_insert_orientation =
                GibsonUiInsertOrientation::from_plan_orientation(&fragment.orientation);
            self.gibson_extra_inserts = plan
                .fragments
                .iter()
                .skip(1)
                .map(|fragment| GibsonUiInsertRow {
                    seq_id: fragment.seq_id.trim().to_string(),
                    orientation: GibsonUiInsertOrientation::from_plan_orientation(
                        &fragment.orientation,
                    ),
                })
                .collect();
        } else {
            self.gibson_insert_seq_id.clear();
            self.gibson_insert_orientation = GibsonUiInsertOrientation::Forward;
        }
        self.gibson_overlap_bp_min = plan
            .validation_policy
            .design_targets
            .overlap_bp_min
            .to_string();
        self.gibson_overlap_bp_max = plan
            .validation_policy
            .design_targets
            .overlap_bp_max
            .to_string();
        self.gibson_minimum_overlap_tm_celsius = format!(
            "{:.1}",
            plan.validation_policy
                .design_targets
                .minimum_overlap_tm_celsius
        );
        self.gibson_priming_tm_min_celsius = format!(
            "{:.1}",
            plan.validation_policy
                .design_targets
                .priming_segment_tm_min_celsius
        );
        self.gibson_priming_tm_max_celsius = format!(
            "{:.1}",
            plan.validation_policy
                .design_targets
                .priming_segment_tm_max_celsius
        );
        self.gibson_priming_length_min_bp = plan
            .validation_policy
            .design_targets
            .priming_segment_min_length_bp
            .to_string();
        self.gibson_priming_length_max_bp = plan
            .validation_policy
            .design_targets
            .priming_segment_max_length_bp
            .to_string();
        self.gibson_unique_restriction_site_enzyme_name = plan
            .validation_policy
            .desired_unique_restriction_site_enzyme_name
            .clone()
            .unwrap_or_default();
        self.gibson_output_id_hint = plan.product.output_id_hint.trim().to_string();
    }

    pub(super) fn reopen_gibson_specialist_from_operation(
        &mut self,
        op_id: &str,
    ) -> std::result::Result<bool, String> {
        let plan_json = {
            let engine = self.engine.read().map_err(|_| {
                "Could not read engine state while reopening Gibson plan".to_string()
            })?;
            engine
                .operation_log()
                .iter()
                .find(|record| record.result.op_id == op_id)
                .and_then(|record| match &record.op {
                    Operation::ApplyGibsonAssemblyPlan { plan_json } => Some(plan_json.clone()),
                    _ => None,
                })
        };
        let Some(plan_json) = plan_json else {
            return Ok(false);
        };
        let plan: GibsonAssemblyPlan = serde_json::from_str(&plan_json).map_err(|err| {
            format!("Could not parse stored Gibson plan JSON from '{op_id}': {err}")
        })?;
        self.load_gibson_plan_into_ui(&plan);
        self.gibson_preview_output = None;
        self.gibson_preview_svg_uri.clear();
        self.gibson_status = format!(
            "Loaded Gibson cloning operation '{}' back into the specialist. Review the design or apply again explicitly.",
            op_id
        );
        self.open_gibson_dialog();
        self.run_gibson_preview();
        Ok(true)
    }

    pub(super) fn run_gibson_apply(&mut self) {
        match self.gibson_preview_output.as_ref() {
            Some(preview) if !preview.can_execute => {
                self.gibson_status = Self::gibson_apply_blocked_status(preview);
                return;
            }
            None => {
                self.gibson_status =
                    "Run Gibson preview first so apply uses one reviewed, executable design."
                        .to_string();
                return;
            }
            Some(_) => {}
        }
        let plan = match self.build_gibson_plan_from_ui() {
            Ok(plan) => plan,
            Err(err) => {
                self.gibson_status = format!("Gibson apply blocked: {err}");
                return;
            }
        };
        let request_json = match serde_json::to_string_pretty(&plan) {
            Ok(text) => text,
            Err(err) => {
                self.gibson_status = format!("Could not serialize Gibson plan JSON: {err}");
                return;
            }
        };
        let command = ShellCommand::GibsonApply { request_json };
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => {
                let result_value = output
                    .get("result")
                    .cloned()
                    .ok_or_else(|| "Missing result payload in Gibson apply output".to_string());
                match result_value.and_then(|value| {
                    serde_json::from_value::<OpResult>(value)
                        .map_err(|err| format!("Could not parse Gibson apply result JSON: {err}"))
                }) {
                    Ok(result) => {
                        let created = if result.created_seq_ids.is_empty() {
                            "-".to_string()
                        } else {
                            result.created_seq_ids.join(", ")
                        };
                        self.gibson_status = format!(
                            "Applied Gibson cloning via {}. Created nodes: {}. Click this Gibson operation in lineage to reopen the specialist later.",
                            result.op_id, created
                        );
                    }
                    Err(err) => {
                        self.gibson_status = err;
                    }
                }
            }
            Err(err) => {
                self.gibson_status = format!("Gibson apply failed: {err}");
            }
        }
    }

    pub(super) fn render_gibson_preview_svg(
        &mut self,
        preview: &GibsonAssemblyPreview,
    ) -> std::result::Result<(), String> {
        let spec = Self::resolve_gibson_preview_cartoon_spec(preview)?;
        let svg = render_protocol_cartoon_spec_svg(&spec);
        let stamp = now_unix_ms_u64();
        let svg_path = env::temp_dir().join(format!("gentle_gibson_preview_{stamp}.svg"));
        fs::write(&svg_path, &svg).map_err(|e| {
            format!(
                "Could not write temporary Gibson preview SVG '{}': {e}",
                svg_path.display()
            )
        })?;

        let mut opt = usvg::Options {
            resources_dir: svg_path.parent().map(|path| path.to_path_buf()),
            ..usvg::Options::default()
        };
        opt.fontdb_mut().load_system_fonts();
        let tree = usvg::Tree::from_str(&svg, &opt)
            .map_err(|e| format!("Could not parse Gibson preview SVG: {e}"))?;
        let pixmap_size = tree
            .size()
            .to_int_size()
            .scale_by(1.0)
            .ok_or_else(|| "Could not determine Gibson preview raster size".to_string())?;
        let mut pixmap = tiny_skia::Pixmap::new(pixmap_size.width(), pixmap_size.height())
            .ok_or_else(|| {
                format!(
                    "Could not allocate Gibson preview raster {}x{}",
                    pixmap_size.width(),
                    pixmap_size.height()
                )
            })?;
        resvg::render(&tree, tiny_skia::Transform::default(), &mut pixmap.as_mut());
        let png_path = env::temp_dir().join(format!("gentle_gibson_preview_{stamp}.png"));
        pixmap.save_png(&png_path).map_err(|e| {
            format!(
                "Could not write temporary Gibson preview PNG '{}': {e}",
                png_path.display()
            )
        })?;
        self.gibson_preview_svg_uri = Self::file_uri_from_path(&png_path);
        Ok(())
    }

    pub(super) fn resolve_gibson_preview_cartoon_spec(
        preview: &GibsonAssemblyPreview,
    ) -> std::result::Result<crate::protocol_cartoon::ProtocolCartoonSpec, String> {
        if let Some(spec) = preview.cartoon.resolved_spec.as_ref() {
            return Ok(spec.clone());
        }
        let protocol =
            ProtocolCartoonKind::parse_id(&preview.cartoon.protocol_id).ok_or_else(|| {
                format!(
                    "Unknown protocol cartoon id '{}' in Gibson preview",
                    preview.cartoon.protocol_id
                )
            })?;
        let template = protocol_cartoon_template_for_kind(&protocol);
        let mut spec =
            resolve_protocol_cartoon_template_with_bindings(&template, &preview.cartoon.bindings)?;
        if !preview.cartoon.title.trim().is_empty() {
            spec.title = preview.cartoon.title.trim().to_string();
        }
        if !preview.cartoon.summary.trim().is_empty() {
            spec.summary = preview.cartoon.summary.trim().to_string();
        }
        Ok(spec)
    }

    pub(super) fn run_gibson_preview(&mut self) {
        if let Err(err) = self.validate_gibson_preview_readiness() {
            self.gibson_status = format!("Gibson preview blocked: {err}");
            self.gibson_preview_output = None;
            self.gibson_preview_svg_uri.clear();
            return;
        }
        let plan = match self.build_gibson_plan_from_ui() {
            Ok(plan) => plan,
            Err(err) => {
                self.gibson_status = format!("Gibson plan error: {err}");
                self.gibson_preview_output = None;
                self.gibson_preview_svg_uri.clear();
                return;
            }
        };
        let request_json = match serde_json::to_string_pretty(&plan) {
            Ok(text) => text,
            Err(err) => {
                self.gibson_status = format!("Could not serialize Gibson plan JSON: {err}");
                self.gibson_preview_output = None;
                self.gibson_preview_svg_uri.clear();
                return;
            }
        };
        let command = ShellCommand::GibsonPreview {
            request_json,
            output_path: None,
        };
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => match serde_json::from_value::<GibsonAssemblyPreview>(output) {
                Ok(preview) => {
                    self.gibson_status = if preview.can_execute {
                        format!(
                            "Gibson preview ready: {} junction(s), {} primer suggestion(s)",
                            preview.resolved_junctions.len(),
                            preview.primer_suggestions.len()
                        )
                    } else {
                        format!(
                            "Gibson preview produced {} blocking error(s)",
                            preview.errors.len()
                        )
                    };
                    self.gibson_preview_output = Some(preview.clone());
                    if let Err(err) = self.render_gibson_preview_svg(&preview) {
                        let base_status = self.gibson_status.clone();
                        self.gibson_status =
                            format!("{} | preview render note: {err}", base_status);
                    }
                }
                Err(err) => {
                    self.gibson_status =
                        format!("Could not parse Gibson preview response JSON: {err}");
                    self.gibson_preview_output = None;
                    self.gibson_preview_svg_uri.clear();
                }
            },
            Err(err) => {
                self.gibson_status = format!("Gibson preview failed: {err}");
                self.gibson_preview_output = None;
                self.gibson_preview_svg_uri.clear();
            }
        }
    }

    pub(super) fn export_gibson_plan_json(&mut self) {
        let plan = match self.build_gibson_plan_from_ui() {
            Ok(plan) => plan,
            Err(err) => {
                self.gibson_status = format!("Gibson export blocked: {err}");
                return;
            }
        };
        let Some(path) = rfd::FileDialog::new()
            .set_file_name("gibson.plan.json")
            .add_filter("JSON", &["json"])
            .save_file()
        else {
            self.gibson_status = "Gibson plan export canceled".to_string();
            return;
        };
        let mut text = match serde_json::to_string_pretty(&plan) {
            Ok(text) => text,
            Err(err) => {
                self.gibson_status = format!("Could not serialize Gibson plan JSON: {err}");
                return;
            }
        };
        text.push('\n');
        match fs::write(&path, text) {
            Ok(()) => {
                self.gibson_status = format!("Exported Gibson plan JSON to '{}'", path.display())
            }
            Err(err) => {
                self.gibson_status = format!(
                    "Could not export Gibson plan JSON '{}': {err}",
                    path.display()
                )
            }
        }
    }

    pub(super) fn export_gibson_preview_json(&mut self) {
        let Some(preview) = self.gibson_preview_output.as_ref() else {
            self.gibson_status = "Run Gibson preview before exporting preview JSON".to_string();
            return;
        };
        let Some(path) = rfd::FileDialog::new()
            .set_file_name("gibson.preview.json")
            .add_filter("JSON", &["json"])
            .save_file()
        else {
            self.gibson_status = "Gibson preview export canceled".to_string();
            return;
        };
        let mut text = match serde_json::to_string_pretty(preview) {
            Ok(text) => text,
            Err(err) => {
                self.gibson_status = format!("Could not serialize Gibson preview JSON: {err}");
                return;
            }
        };
        text.push('\n');
        match fs::write(&path, text) {
            Ok(()) => {
                self.gibson_status = format!("Exported Gibson preview JSON to '{}'", path.display())
            }
            Err(err) => {
                self.gibson_status = format!(
                    "Could not export Gibson preview JSON '{}': {err}",
                    path.display()
                )
            }
        }
    }

    pub(super) fn export_gibson_primer_summary(&mut self) {
        let Some(preview) = self.gibson_preview_output.as_ref() else {
            self.gibson_status = "Run Gibson preview before exporting primer summary".to_string();
            return;
        };
        let Some(path) = rfd::FileDialog::new()
            .set_file_name("gibson.primers.txt")
            .add_filter("Text", &["txt"])
            .save_file()
        else {
            self.gibson_status = "Gibson primer export canceled".to_string();
            return;
        };
        let mut text = String::new();
        for primer in &preview.primer_suggestions {
            text.push_str(&format!("{}\n", primer.primer_id));
            text.push_str(&format!("  full: {}\n", primer.full_sequence));
            text.push_str(&format!(
                "  5' overlap: {} ({} bp, {:.1} C)\n",
                primer.overlap_5prime.sequence,
                primer.overlap_5prime.length_bp,
                primer.overlap_5prime.tm_celsius
            ));
            text.push_str(&format!(
                "  3' priming: {} ({} bp, {:.1} C, hits={})\n\n",
                primer.priming_3prime.sequence,
                primer.priming_3prime.length_bp,
                primer.priming_3prime.tm_celsius,
                primer.priming_3prime.anneal_hits
            ));
        }
        match fs::write(&path, text) {
            Ok(()) => {
                self.gibson_status =
                    format!("Exported Gibson primer summary to '{}'", path.display())
            }
            Err(err) => {
                self.gibson_status = format!(
                    "Could not export Gibson primer summary '{}': {err}",
                    path.display()
                )
            }
        }
    }

    pub(super) fn export_gibson_cartoon_svg(&mut self) {
        let Some(preview) = self.gibson_preview_output.as_ref() else {
            self.gibson_status = "Run Gibson preview before exporting cartoon SVG".to_string();
            return;
        };
        let Some(path) = rfd::FileDialog::new()
            .set_file_name("gibson.preview.svg")
            .add_filter("SVG", &["svg"])
            .save_file()
        else {
            self.gibson_status = "Gibson cartoon export canceled".to_string();
            return;
        };
        let spec = match Self::resolve_gibson_preview_cartoon_spec(preview) {
            Ok(spec) => spec,
            Err(err) => {
                self.gibson_status = format!("Could not resolve Gibson cartoon bindings: {err}");
                return;
            }
        };
        let svg = render_protocol_cartoon_spec_svg(&spec);
        match fs::write(&path, svg) {
            Ok(()) => {
                self.gibson_status = format!("Exported Gibson cartoon SVG to '{}'", path.display())
            }
            Err(err) => {
                self.gibson_status = format!(
                    "Could not export Gibson cartoon SVG '{}': {err}",
                    path.display()
                )
            }
        }
    }

    pub(super) fn handoff_gibson_preview_to_routine_assistant(&mut self) {
        let Some(preview) = self.gibson_preview_output.clone() else {
            self.gibson_status =
                "Run Gibson preview before opening Routine Assistant handoff".to_string();
            return;
        };
        self.open_routine_assistant_dialog();
        self.refresh_routine_assistant_candidates();
        self.routine_assistant_goal = "Gibson assembly".to_string();
        self.routine_assistant_query = "gibson".to_string();
        self.routine_assistant_selected_routine_id = preview.routine_handoff.routine_id.clone();
        self.routine_assistant_compare_routine_id.clear();
        self.routine_assistant_bindings = preview.routine_handoff.bindings.clone();
        self.routine_assistant_explain_output = None;
        self.routine_assistant_compare_output = None;
        self.routine_assistant_preflight_output = None;
        self.routine_assistant_execute_output = None;
        if preview.routine_handoff.supported {
            self.routine_assistant_stage = RoutineAssistantStage::Parameters;
            self.routine_assistant_status = format!(
                "Gibson preview handed off to '{}' with {} bound input(s)",
                preview.routine_handoff.routine_id,
                preview.routine_handoff.bindings.len()
            );
        } else {
            self.routine_assistant_stage = RoutineAssistantStage::GoalAndCandidates;
            self.routine_assistant_status = preview.routine_handoff.reason.unwrap_or_else(|| {
                "Gibson preview is currently preview-only for this opening mode".to_string()
            });
        }
    }

    pub(super) fn current_gibson_destination_opening_suggestions(
        &self,
    ) -> std::result::Result<Vec<GibsonDestinationOpeningSuggestion>, String> {
        let destination_id = self.gibson_destination_seq_id.trim();
        if destination_id.is_empty() {
            return Ok(vec![]);
        }
        let engine = self.engine.read().map_err(|_| {
            "Could not read engine state for Gibson opening suggestions".to_string()
        })?;
        suggest_gibson_destination_openings(&engine, destination_id).map_err(|err| err.message)
    }

    pub(super) fn gibson_destination_suggestion_geometry_label(
        suggestion: &GibsonDestinationOpeningSuggestion,
    ) -> String {
        match suggestion.end_geometry.as_str() {
            "blunt" => "blunt".to_string(),
            "5prime_overhang" => {
                format!("5' overhang ({} bp)", suggestion.overhang_bp.unwrap_or(0))
            }
            "3prime_overhang" => {
                format!("3' overhang ({} bp)", suggestion.overhang_bp.unwrap_or(0))
            }
            "feature_span" => "region".to_string(),
            _ => suggestion.end_geometry.clone(),
        }
    }

    pub(super) fn humanize_gibson_ui_text(raw: &str) -> String {
        raw.replace(" C", " °C")
    }

    pub(super) fn gibson_tm_label(
        ui: &mut Ui,
        prefix: &str,
        suffix: &str,
        text_style: egui::TextStyle,
        strong: bool,
    ) -> egui::Response {
        let font_id = text_style.resolve(ui.style());
        let color = if strong {
            ui.visuals().strong_text_color()
        } else {
            ui.visuals().text_color()
        };
        let mut job = egui::text::LayoutJob::default();
        let base_format = egui::TextFormat {
            font_id: font_id.clone(),
            color,
            valign: egui::Align::Center,
            ..Default::default()
        };
        let subscript_format = egui::TextFormat {
            font_id: egui::FontId::new(font_id.size * 0.72, font_id.family.clone()),
            color,
            valign: egui::Align::BOTTOM,
            ..Default::default()
        };
        if !prefix.is_empty() {
            job.append(prefix, 0.0, base_format.clone());
        }
        job.append("T", 0.0, base_format.clone());
        job.append("m", 0.0, subscript_format);
        if !suffix.is_empty() {
            job.append(suffix, 0.0, base_format);
        }
        ui.label(job)
    }

    pub(super) fn gibson_rich_text_job(
        ui: &Ui,
        raw: &str,
        text_style: egui::TextStyle,
        color: egui::Color32,
    ) -> egui::text::LayoutJob {
        let text = Self::humanize_gibson_ui_text(raw);
        let font_id = text_style.resolve(ui.style());
        let base_format = egui::TextFormat {
            font_id: font_id.clone(),
            color,
            valign: egui::Align::Center,
            ..Default::default()
        };
        let subscript_format = egui::TextFormat {
            font_id: egui::FontId::new(font_id.size * 0.72, font_id.family.clone()),
            color,
            valign: egui::Align::BOTTOM,
            ..Default::default()
        };
        let mut job = egui::text::LayoutJob::default();
        let mut remaining = text.as_str();
        while let Some(idx) = remaining.find("Tm") {
            let (before, after_before) = remaining.split_at(idx);
            if !before.is_empty() {
                job.append(before, 0.0, base_format.clone());
            }
            job.append("T", 0.0, base_format.clone());
            job.append("m", 0.0, subscript_format.clone());
            remaining = &after_before[2..];
        }
        if !remaining.is_empty() {
            job.append(remaining, 0.0, base_format);
        }
        job
    }

    pub(super) fn gibson_rich_text_label(
        ui: &mut Ui,
        raw: &str,
        text_style: egui::TextStyle,
        color: egui::Color32,
    ) -> egui::Response {
        ui.label(Self::gibson_rich_text_job(ui, raw, text_style, color))
    }

    pub(super) fn gibson_active_context_summary(&self) -> String {
        match self.active_dna_window_context() {
            Some((seq_id, Some((start, end)))) => {
                format!("Active DNA window: {seq_id} | selection {start}..{end} (0-based)")
            }
            Some((seq_id, None)) => format!("Active DNA window: {seq_id} | no selection"),
            None => "Active DNA window: none".to_string(),
        }
    }

    pub(super) fn gibson_use_active_sequence_tooltip(&self) -> String {
        match self.active_dna_window_context() {
            Some((seq_id, _)) => format!(
                "Fill the destination from the currently focused DNA window.\nCurrent active sequence: {seq_id}"
            ),
            None => "Fill the destination from the currently focused DNA window.\nNo active DNA window is currently available.".to_string(),
        }
    }

    pub(super) fn gibson_use_active_selection_tooltip(&self) -> String {
        match self.active_dna_window_context() {
            Some((seq_id, Some((start, end)))) => format!(
                "Fill the Gibson cut/opening edges from the current selection in the focused DNA window.\nCurrent active selection: {seq_id} {start}..{end} (0-based)"
            ),
            Some((seq_id, None)) => format!(
                "Fill the Gibson cut/opening edges from the current selection in the focused DNA window.\nCurrent active DNA window: {seq_id} (no active selection)"
            ),
            None => "Fill the Gibson cut/opening edges from the current selection in the focused DNA window.\nNo active DNA window is currently available.".to_string(),
        }
    }

    pub(super) fn gibson_destination_suggestion_help_text(
        suggestion: &GibsonDestinationOpeningSuggestion,
    ) -> String {
        let why = if suggestion.in_mcs_context {
            "This unique cutter is named in the selected MCS annotation, so it is shown first."
        } else {
            "This unique cutter is available on the destination, but it is not part of the selected MCS annotation."
        };
        format!("{}\n{}", suggestion.summary, why)
    }

    pub(super) fn gibson_destination_suggestion_display_label(
        suggestion: &GibsonDestinationOpeningSuggestion,
    ) -> String {
        if suggestion.in_mcs_context {
            format!("{} (MCS)", suggestion.label)
        } else {
            suggestion.label.clone()
        }
    }

    pub(super) fn gibson_destination_suggestion_cut_label(
        suggestion: &GibsonDestinationOpeningSuggestion,
    ) -> String {
        suggestion.rebase_cut_summary.clone()
    }

    pub(super) fn gibson_destination_suggestion_feature_label(
        suggestion: &GibsonDestinationOpeningSuggestion,
    ) -> String {
        if suggestion.feature_context.trim().is_empty() {
            "-".to_string()
        } else {
            suggestion.feature_context.clone()
        }
    }

    pub(super) fn open_gibson_show_all_unique_cutters(&mut self) {
        self.gibson_show_all_unique_cutters = true;
        self.gibson_status =
            "Searching all known unique cleavage sites on the current destination".to_string();
    }

    pub(super) fn compact_gibson_sequence(sequence: &str, edge_bp: usize) -> String {
        if sequence.len() <= edge_bp.saturating_mul(2).saturating_add(3) {
            sequence.to_string()
        } else {
            format!(
                "{}...{}",
                &sequence[..edge_bp],
                &sequence[sequence.len().saturating_sub(edge_bp)..]
            )
        }
    }

    pub(super) fn current_gibson_destination_sequence_record(&self) -> Option<DNAsequence> {
        let destination_id = self.gibson_destination_seq_id.trim().to_string();
        if destination_id.is_empty() {
            return None;
        }
        self.engine
            .read()
            .ok()
            .and_then(|engine| engine.state().sequences.get(&destination_id).cloned())
    }

    pub(super) fn validate_gibson_preview_readiness(&self) -> std::result::Result<(), String> {
        if self.gibson_destination_seq_id.trim().is_empty() {
            return Err("Choose a destination sequence first.".to_string());
        }
        if self.gibson_ui_insert_rows().is_empty() {
            return Err("Choose at least one insert sequence first.".to_string());
        }
        if self.gibson_opening_mode == GibsonUiOpeningMode::DefinedSite {
            if self.gibson_opening_start_0based.trim().is_empty()
                || self.gibson_opening_end_0based_exclusive.trim().is_empty()
            {
                return Err(
                    "Set both left and right cut edges or choose a cutter suggestion first."
                        .to_string(),
                );
            }
            if self
                .gibson_opening_start_0based
                .trim()
                .parse::<usize>()
                .is_err()
            {
                return Err("Left cut edge must be a non-negative integer.".to_string());
            }
            if self
                .gibson_opening_end_0based_exclusive
                .trim()
                .parse::<usize>()
                .is_err()
            {
                return Err("Right cut edge must be a non-negative integer.".to_string());
            }
        }
        Ok(())
    }

    pub(super) fn current_gibson_opening_edges(&self) -> Option<(usize, usize)> {
        let start = self
            .gibson_opening_start_0based
            .trim()
            .parse::<usize>()
            .ok()?;
        let end = self
            .gibson_opening_end_0based_exclusive
            .trim()
            .parse::<usize>()
            .ok()?;
        Some((start, end))
    }

    pub(super) fn current_gibson_opening_suggestion_for_edges<'a>(
        suggestions: &'a [GibsonDestinationOpeningSuggestion],
        start: usize,
        end: usize,
    ) -> Option<&'a GibsonDestinationOpeningSuggestion> {
        suggestions
            .iter()
            .find(|row| row.start_0based == start && row.end_0based_exclusive == end)
    }

    pub(super) fn gibson_linear_suffix(
        sequence: &str,
        end_exclusive: usize,
        len: usize,
    ) -> Option<String> {
        if end_exclusive > sequence.len() {
            return None;
        }
        let start = end_exclusive.saturating_sub(len);
        sequence.get(start..end_exclusive).map(str::to_string)
    }

    pub(super) fn gibson_linear_prefix(sequence: &str, start: usize, len: usize) -> Option<String> {
        if start > sequence.len() {
            return None;
        }
        let end = (start + len).min(sequence.len());
        sequence.get(start..end).map(str::to_string)
    }

    pub(super) fn gibson_circular_slice(
        sequence: &str,
        start: usize,
        len: usize,
    ) -> Option<String> {
        if sequence.is_empty() || len == 0 || start > sequence.len() {
            return None;
        }
        let total = sequence.len();
        let start = start % total;
        let mut out = String::with_capacity(len);
        for idx in 0..len {
            let pos = (start + idx) % total;
            out.push(*sequence.as_bytes().get(pos)? as char);
        }
        Some(out)
    }

    pub(super) fn gibson_suffix_before_edge(
        sequence: &str,
        edge: usize,
        len: usize,
        circular: bool,
    ) -> Option<String> {
        if len == 0 {
            return Some(String::new());
        }
        if circular {
            if sequence.is_empty() || edge > sequence.len() {
                return None;
            }
            let total = sequence.len();
            let start = (edge + total - (len % total)) % total;
            Self::gibson_circular_slice(sequence, start, len)
        } else {
            Self::gibson_linear_suffix(sequence, edge, len)
        }
    }

    pub(super) fn gibson_prefix_after_edge(
        sequence: &str,
        edge: usize,
        len: usize,
        circular: bool,
    ) -> Option<String> {
        if len == 0 {
            return Some(String::new());
        }
        if circular {
            Self::gibson_circular_slice(sequence, edge, len)
        } else {
            Self::gibson_linear_prefix(sequence, edge, len)
        }
    }

    pub(super) fn gibson_between_edges(
        sequence: &str,
        start: usize,
        end: usize,
        circular: bool,
    ) -> Option<String> {
        if start > sequence.len() || end > sequence.len() {
            return None;
        }
        if start == end {
            return Some(String::new());
        }
        if circular && start > end {
            let mut out = String::new();
            out.push_str(sequence.get(start..sequence.len())?);
            out.push_str(sequence.get(0..end)?);
            return Some(out);
        }
        sequence.get(start..end).map(str::to_string)
    }

    pub(super) fn gibson_opening_detail_text(
        &self,
        suggestions: &[GibsonDestinationOpeningSuggestion],
    ) -> Option<String> {
        let destination = self.current_gibson_destination_sequence_record()?;
        let sequence = destination.get_forward_string().to_ascii_uppercase();
        let (start, end) = self.current_gibson_opening_edges()?;
        if start > sequence.len() || end > sequence.len() {
            return None;
        }
        let matched = Self::current_gibson_opening_suggestion_for_edges(suggestions, start, end);
        let mut lines = vec![];
        if let Some(suggestion) = matched {
            lines.push(format!(
                "{}: {}",
                Self::gibson_destination_suggestion_display_label(suggestion),
                suggestion.rebase_cut_summary
            ));
        }
        let context_bp = 12;
        let left = Self::gibson_suffix_before_edge(
            &sequence,
            start,
            context_bp,
            destination.is_circular(),
        )?;
        let middle = Self::gibson_between_edges(&sequence, start, end, destination.is_circular())?;
        let right =
            Self::gibson_prefix_after_edge(&sequence, end, context_bp, destination.is_circular())?;
        let left_ellipsis = if destination.is_circular() || start > context_bp {
            "..."
        } else {
            ""
        };
        let right_ellipsis = if destination.is_circular() || end + context_bp < sequence.len() {
            "..."
        } else {
            ""
        };
        if middle.is_empty() {
            lines.push(format!(
                "destination cut: {left_ellipsis}{left}|{right}{right_ellipsis}"
            ));
        } else {
            lines.push(format!(
                "destination opening: {left_ellipsis}{left}|{middle}|{right}{right_ellipsis}"
            ));
        }
        Some(lines.join("\n"))
    }

    pub(super) fn gibson_opened_destination_arm_text(&self) -> Option<String> {
        let destination = self.current_gibson_destination_sequence_record()?;
        let sequence = destination.get_forward_string().to_ascii_uppercase();
        let (start, end) = self.current_gibson_opening_edges()?;
        let context_bp = 12;
        let left_context = Self::gibson_suffix_before_edge(
            &sequence,
            start,
            context_bp,
            destination.is_circular(),
        )?;
        let right_context =
            Self::gibson_prefix_after_edge(&sequence, end, context_bp, destination.is_circular())?;
        Some(format!(
            "left arm  : ...{}|\nright arm : |{}...",
            left_context, right_context
        ))
    }

    pub(super) fn gibson_resolved_destination_arm_text(
        &self,
        preview: &GibsonAssemblyPreview,
    ) -> Option<String> {
        let destination = self.current_gibson_destination_sequence_record()?;
        let sequence = destination.get_forward_string().to_ascii_uppercase();
        let start = preview.destination.opening_start_0based?;
        let end = preview.destination.opening_end_0based_exclusive?;
        let left_junction = preview
            .resolved_junctions
            .iter()
            .find(|junction| junction.left_member_id == preview.destination.left_end_id)?;
        let right_junction = preview
            .resolved_junctions
            .iter()
            .find(|junction| junction.right_member_id == preview.destination.right_end_id)?;
        let context_bp = 8;
        let left_context = Self::gibson_suffix_before_edge(
            &sequence,
            start,
            context_bp,
            destination.is_circular(),
        )?;
        let right_context =
            Self::gibson_prefix_after_edge(&sequence, end, context_bp, destination.is_circular())?;
        Some(format!(
            "left 5' arm ({} bp overlap): ...{}[{}]|\nright 5' arm ({} bp overlap): |[{}]{}...",
            left_junction.overlap_bp,
            left_context,
            Self::compact_gibson_sequence(&left_junction.overlap_sequence, 8),
            right_junction.overlap_bp,
            Self::compact_gibson_sequence(&right_junction.overlap_sequence, 8),
            right_context
        ))
    }

    pub(super) fn gibson_primer_construction_text(
        preview: &GibsonAssemblyPreview,
    ) -> Option<String> {
        if preview.primer_suggestions.is_empty() {
            return None;
        }
        Some(
            preview
                .primer_suggestions
                .iter()
                .map(|primer| {
                    format!(
                        "{}\n  full       : {}\n  5' overlap : {}\n  3' priming : {}",
                        primer.primer_id,
                        Self::compact_gibson_sequence(&primer.full_sequence, 10),
                        Self::compact_gibson_sequence(&primer.overlap_5prime.sequence, 8),
                        Self::compact_gibson_sequence(&primer.priming_3prime.sequence, 8),
                    )
                })
                .collect::<Vec<_>>()
                .join("\n\n"),
        )
    }

    pub(super) fn displayed_gibson_destination_opening_suggestions(
        suggestions: &[GibsonDestinationOpeningSuggestion],
        show_all_unique_cutters: bool,
    ) -> (Vec<GibsonDestinationOpeningSuggestion>, usize) {
        let has_mcs_linked = suggestions.iter().any(|row| row.in_mcs_context);
        if show_all_unique_cutters || !has_mcs_linked {
            return (suggestions.to_vec(), 0);
        }
        let visible = suggestions
            .iter()
            .filter(|row| row.in_mcs_context)
            .cloned()
            .collect::<Vec<_>>();
        let hidden_count = suggestions.len().saturating_sub(visible.len());
        (visible, hidden_count)
    }

    pub(super) fn gibson_preview_findings_text(preview: &GibsonAssemblyPreview) -> String {
        let mut lines = vec![];
        if !preview.errors.is_empty() {
            lines.push(format!("blocking errors: {}", preview.errors.len()));
            for row in &preview.errors {
                lines.push(format!("- {}", Self::humanize_gibson_ui_text(row)));
            }
        }
        if !preview.warnings.is_empty() {
            lines.push(format!("warnings: {}", preview.warnings.len()));
            for row in &preview.warnings {
                lines.push(format!("- {}", Self::humanize_gibson_ui_text(row)));
            }
        }
        if !preview.notes.is_empty() {
            lines.push("notes:".to_string());
            for row in &preview.notes {
                lines.push(format!("- {}", Self::humanize_gibson_ui_text(row)));
            }
        }
        lines.join("\n")
    }

    pub(super) fn is_gibson_priming_error(row: &str) -> bool {
        row.contains("Could not derive a Gibson priming segment")
    }

    pub(super) fn is_gibson_overlap_error(row: &str) -> bool {
        row.contains("Could not derive an overlap length")
            || row.contains("overlap Tm")
            || row.contains("Terminal overlap regions")
            || row.contains("requires overlap_partition")
            || row.contains("does not match required_overlap_bp")
            || row.contains("not destination-derived")
            || row.contains("Destination does not have enough flank sequence")
    }

    pub(super) fn gibson_target_review_rows(
        &self,
        preview: &GibsonAssemblyPreview,
    ) -> Vec<(egui::Color32, String)> {
        let ok_color = egui::Color32::from_rgb(60, 130, 75);
        let warn_color = egui::Color32::from_rgb(170, 120, 20);
        let error_color = egui::Color32::from_rgb(190, 70, 70);
        let overlap_error_count = preview
            .errors
            .iter()
            .filter(|row| Self::is_gibson_overlap_error(row))
            .count();
        let priming_error_count = preview
            .errors
            .iter()
            .filter(|row| Self::is_gibson_priming_error(row))
            .count();
        let expected_junction_count = preview.fragments.len().saturating_add(1);
        let primed_fragment_count = preview
            .primer_suggestions
            .iter()
            .map(|row| row.fragment_id.clone())
            .collect::<HashSet<_>>()
            .len();

        let overlap_line = if preview.resolved_junctions.len() == expected_junction_count
            && overlap_error_count == 0
        {
            (
                ok_color,
                format!(
                    "Overlap target: resolved {}/{} junction overlaps cleanly.",
                    preview.resolved_junctions.len(),
                    expected_junction_count
                ),
            )
        } else {
            (
                error_color,
                format!(
                    "Overlap target: only {}/{} junction overlaps resolve under the current settings.",
                    preview.resolved_junctions.len(),
                    expected_junction_count
                ),
            )
        };

        let priming_line = if primed_fragment_count == preview.fragments.len()
            && priming_error_count == 0
        {
            (
                ok_color,
                format!(
                    "PCR priming target: complete primer pairs were derived for all {}/{} insert fragments.",
                    primed_fragment_count,
                    preview.fragments.len()
                ),
            )
        } else if preview.resolved_junctions.len() == expected_junction_count
            && overlap_error_count == 0
        {
            (
                warn_color,
                format!(
                    "PCR priming target: complete primer pairs were derived for {}/{} insert fragments. The current blockers are in the 3' priming window, not in the 5' Gibson overlaps.",
                    primed_fragment_count,
                    preview.fragments.len()
                ),
            )
        } else {
            (
                error_color,
                format!(
                    "PCR priming target: complete primer pairs were derived for {}/{} insert fragments.",
                    primed_fragment_count,
                    preview.fragments.len()
                ),
            )
        };

        let mut rows = vec![overlap_line, priming_line];
        if priming_error_count > 0
            && preview.resolved_junctions.len() == expected_junction_count
            && overlap_error_count == 0
        {
            let adjustment_text = if preview.suggested_design_adjustments.is_empty() {
                format!(
                    "try increasing max priming length above {} bp or lowering the minimum priming Tm below {} °C",
                    self.gibson_priming_length_max_bp.trim(),
                    self.gibson_priming_tm_min_celsius.trim()
                )
            } else {
                preview
                    .suggested_design_adjustments
                    .iter()
                    .map(Self::gibson_design_adjustment_brief)
                    .collect::<Vec<_>>()
                    .join(" or ")
            };
            rows.push((
                warn_color,
                format!("Design hint: if biologically acceptable, {adjustment_text}."),
            ));
        }
        rows
    }

    pub(super) fn gibson_design_adjustment_brief(
        adjustment: &GibsonSuggestedDesignAdjustment,
    ) -> String {
        match adjustment.target {
            GibsonDesignAdjustmentTarget::PrimingSegmentMaxLengthBp => format!(
                "set max priming length to {:.0} bp",
                adjustment.suggested_value.round()
            ),
            GibsonDesignAdjustmentTarget::PrimingSegmentTmMinCelsius => {
                format!("set min priming Tm to {:.1} °C", adjustment.suggested_value)
            }
        }
    }

    pub(super) fn apply_gibson_design_adjustment(
        &mut self,
        adjustment: &GibsonSuggestedDesignAdjustment,
    ) {
        match adjustment.target {
            GibsonDesignAdjustmentTarget::PrimingSegmentMaxLengthBp => {
                self.gibson_priming_length_max_bp =
                    format!("{:.0}", adjustment.suggested_value.round());
            }
            GibsonDesignAdjustmentTarget::PrimingSegmentTmMinCelsius => {
                self.gibson_priming_tm_min_celsius = format!("{:.1}", adjustment.suggested_value);
            }
        }
        self.run_gibson_preview();
    }

    pub(super) fn gibson_preview_junctions_text(preview: &GibsonAssemblyPreview) -> String {
        preview
            .resolved_junctions
            .iter()
            .map(|junction| {
                format!(
                    "{}\n  members: {} -> {}\n  overlap: {} bp | {}\n  Tm: {:.1} °C\n  source: {}",
                    junction.junction_id,
                    junction.left_member_id,
                    junction.right_member_id,
                    junction.overlap_bp,
                    junction.overlap_sequence,
                    junction.overlap_tm_celsius,
                    junction.overlap_source
                )
            })
            .collect::<Vec<_>>()
            .join("\n\n")
    }

    pub(super) fn gibson_preview_primers_text(preview: &GibsonAssemblyPreview) -> String {
        preview
            .primer_suggestions
            .iter()
            .map(|primer| {
                format!(
                    "{}\n  full: {}\n  5' overlap: {} ({} bp, {:.1} °C)\n  3' priming: {} ({} bp, {:.1} °C, hits={})",
                    primer.primer_id,
                    primer.full_sequence,
                    primer.overlap_5prime.sequence,
                    primer.overlap_5prime.length_bp,
                    primer.overlap_5prime.tm_celsius,
                    primer.priming_3prime.sequence,
                    primer.priming_3prime.length_bp,
                    primer.priming_3prime.tm_celsius,
                    primer.priming_3prime.anneal_hits
                )
            })
            .collect::<Vec<_>>()
            .join("\n\n")
    }

    pub(super) fn render_gibson_preview_findings_rich(
        ui: &mut Ui,
        preview: &GibsonAssemblyPreview,
    ) {
        let error_color = egui::Color32::from_rgb(190, 70, 70);
        let warning_color = egui::Color32::from_rgb(170, 120, 20);
        let note_color = ui.visuals().text_color();
        if !preview.errors.is_empty() {
            Self::gibson_rich_text_label(
                ui,
                &format!("blocking errors: {}", preview.errors.len()),
                egui::TextStyle::Small,
                error_color,
            );
            for row in &preview.errors {
                ui.horizontal_wrapped(|ui| {
                    ui.colored_label(error_color, "-");
                    Self::gibson_rich_text_label(ui, row, egui::TextStyle::Small, error_color);
                });
            }
        }
        if !preview.warnings.is_empty() {
            Self::gibson_rich_text_label(
                ui,
                &format!("warnings: {}", preview.warnings.len()),
                egui::TextStyle::Small,
                warning_color,
            );
            for row in &preview.warnings {
                ui.horizontal_wrapped(|ui| {
                    ui.colored_label(warning_color, "-");
                    Self::gibson_rich_text_label(ui, row, egui::TextStyle::Small, warning_color);
                });
            }
        }
        if !preview.notes.is_empty() {
            Self::gibson_rich_text_label(ui, "notes:", egui::TextStyle::Small, note_color);
            for row in &preview.notes {
                ui.horizontal_wrapped(|ui| {
                    ui.colored_label(note_color, "-");
                    Self::gibson_rich_text_label(ui, row, egui::TextStyle::Small, note_color);
                });
            }
        }
    }

    pub(super) fn apply_gibson_destination_opening_suggestion(
        &mut self,
        suggestion: &GibsonDestinationOpeningSuggestion,
    ) {
        self.gibson_opening_mode = GibsonUiOpeningMode::DefinedSite;
        self.gibson_opening_start_0based = suggestion.start_0based.to_string();
        self.gibson_opening_end_0based_exclusive = suggestion.end_0based_exclusive.to_string();
        self.gibson_status = if suggestion.start_0based == suggestion.end_0based_exclusive {
            format!(
                "Applied Gibson opening suggestion: {} -> blunt cutpoint at edge {}",
                suggestion.label, suggestion.start_0based
            )
        } else {
            format!(
                "Applied Gibson opening suggestion: {} -> cut/opening edges {}..{}",
                suggestion.label, suggestion.start_0based, suggestion.end_0based_exclusive
            )
        };
    }

    pub(super) fn cancel_gibson_dialog(&mut self) {
        self.show_gibson_dialog = false;
    }

    pub(super) fn is_lineage_operation_hub(row: &LineageRow) -> bool {
        row.origin == "OperationHub"
    }

    pub(super) fn is_lineage_operation_hub_connector(op_ids: &[String]) -> bool {
        !op_ids.is_empty()
            && op_ids
                .iter()
                .all(|op_id| op_id.ends_with("::hub_in") || op_id.ends_with("::hub_out"))
    }

    pub(super) fn render_gibson_contents(&mut self, ui: &mut Ui) -> bool {
        let mut cancel_clicked = false;
        let close_hover = Self::specialist_window_close_hover_text("Gibson");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            cancel_clicked = true;
        }
        ui.label(
            "Destination-first Gibson specialist: choose a destination, define the opening, choose one or more ordered inserts, then derive overlaps, primers, validation findings, and the factual protocol cartoon from one shared preview path.",
        );
        ui.small(
            "This specialist stays Gibson-specific on purpose: overlap/Tm targets are exposed directly, while generic PCR/qPCR controls remain elsewhere.",
        );
        ui.small(format!(
            "Current limitation: {}",
            Self::gibson_multi_insert_defined_opening_note()
        ));

        let sequence_ids = self.project_sequence_ids_for_blast();
        let active_context_summary = self.gibson_active_context_summary();
        let active_sequence_tooltip = self.gibson_use_active_sequence_tooltip();
        let active_selection_tooltip = self.gibson_use_active_selection_tooltip();
        let destination_topology = {
            let destination_id = self.gibson_destination_seq_id.trim().to_string();
            self.engine
                .read()
                .ok()
                .and_then(|engine| engine.state().sequences.get(&destination_id).cloned())
                .map(|dna| {
                    format!(
                        "{} | {} bp",
                        if dna.is_circular() {
                            "circular"
                        } else {
                            "linear"
                        },
                        dna.len()
                    )
                })
                .unwrap_or_else(|| "not found".to_string())
        };
        let destination_opening_suggestions = self.current_gibson_destination_opening_suggestions();
        let mut selected_opening_suggestion: Option<GibsonDestinationOpeningSuggestion> = None;
        let tm_help = GentleEngine::primer_tm_model_description();
        let preview_readiness = self.validate_gibson_preview_readiness();
        let preview_ready = preview_readiness.is_ok();
        let preview_ready_reason = preview_readiness.err();

        ui.separator();
        ui.heading("Destination");
        ui.horizontal(|ui| {
            ui.label("sequence");
            egui::ComboBox::from_id_salt("gibson_destination_seq_id")
                .selected_text(if self.gibson_destination_seq_id.trim().is_empty() {
                    "(choose destination)"
                } else {
                    self.gibson_destination_seq_id.as_str()
                })
                .show_ui(ui, |ui| {
                    for seq_id in &sequence_ids {
                        ui.selectable_value(
                            &mut self.gibson_destination_seq_id,
                            seq_id.clone(),
                            seq_id,
                        );
                    }
                });
            ui.text_edit_singleline(&mut self.gibson_destination_seq_id);
            if ui
                .button("Use Active Sequence")
                .on_hover_text(active_sequence_tooltip)
                .clicked()
            {
                if let Some((seq_id, _)) = self.active_dna_window_context() {
                    self.gibson_destination_seq_id = seq_id;
                } else {
                    self.gibson_status =
                        "No active DNA window found for Gibson destination prefill".to_string();
                }
            }
        });
        ui.small(format!("topology: {destination_topology}"));
        ui.small(active_context_summary.clone());
        ui.horizontal(|ui| {
            ui.label("opening");
            egui::ComboBox::from_id_salt("gibson_opening_mode")
                .selected_text(self.gibson_opening_mode.label())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.gibson_opening_mode,
                        GibsonUiOpeningMode::DefinedSite,
                        GibsonUiOpeningMode::DefinedSite.label(),
                    );
                    ui.selectable_value(
                        &mut self.gibson_opening_mode,
                        GibsonUiOpeningMode::ExistingTermini,
                        GibsonUiOpeningMode::ExistingTermini.label(),
                    );
                });
            if ui
                .button("Use Active Selection")
                .on_hover_text(active_selection_tooltip)
                .clicked()
            {
                if let Some((seq_id, Some((start, end)))) = self.active_dna_window_context() {
                    self.gibson_destination_seq_id = seq_id;
                    self.gibson_opening_mode = GibsonUiOpeningMode::DefinedSite;
                    self.gibson_opening_start_0based = start.to_string();
                    self.gibson_opening_end_0based_exclusive = end.to_string();
                } else {
                    self.gibson_status =
                        "No active DNA selection found for Gibson opening prefill".to_string();
                }
            }
        });
        ui.horizontal(|ui| {
            ui.label("left cut edge").on_hover_text(
                "0-based edge coordinate between bases on the left side of the Gibson opening.",
            );
            ui.add(
                egui::TextEdit::singleline(&mut self.gibson_opening_start_0based)
                    .desired_width(100.0),
            );
            ui.label("right cut edge").on_hover_text(
                "0-based edge coordinate between bases on the right side of the Gibson opening.",
            );
            ui.add(
                egui::TextEdit::singleline(&mut self.gibson_opening_end_0based_exclusive)
                    .desired_width(120.0),
            );
        });
        ui.small(
            "Equal left/right edges mean one blunt cutpoint. Different left/right edges mean the two primer-relevant termini of a sticky cut; Gibson chew-back removes the opposing 3' ends.",
        );
        match &destination_opening_suggestions {
            Ok(suggestions) if !suggestions.is_empty() => {
                let (displayed_suggestions, hidden_unique_cutters) =
                    Self::displayed_gibson_destination_opening_suggestions(
                        suggestions,
                        self.gibson_show_all_unique_cutters,
                    );
                ui.small(
                    "Suggested openings are unique restriction cutters. Enzymes named by the selected MCS annotation are shown first; other unique cutters can be revealed on demand.",
                );
                let table_height = ((displayed_suggestions.len() + 1) as f32
                    * ui.spacing().interact_size.y)
                    + 10.0;
                const SUGGESTION_COL_WIDTH: f32 = 120.0;
                const ENDS_COL_WIDTH: f32 = 116.0;
                const FEATURE_COL_WIDTH: f32 = 176.0;
                const CUT_COL_WIDTH: f32 = 150.0;
                egui::ScrollArea::horizontal()
                    .id_salt("gibson_destination_opening_suggestions_scroll")
                    .auto_shrink([false, false])
                    .max_height(table_height)
                    .show(ui, |ui| {
                        egui::Grid::new("gibson_destination_opening_suggestions")
                            .num_columns(5)
                            .spacing([8.0, 4.0])
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Suggestion").on_hover_text(
                                    "Unique restriction cutters on the destination sequence. `(MCS)` marks enzymes named by the selected MCS annotation.",
                                );
                                ui.strong("Ends").on_hover_text(
                                    "Blunt / 5' / 3' rows come from actual restriction-cutter geometry.",
                                );
                                ui.strong("Feature").on_hover_text(
                                    "Gene name when available, otherwise another overlapping feature label. MCS-linked cutters also show the MCS location here.",
                                );
                                ui.strong("Cut").on_hover_text(
                                    "Compact REBASE view: motif plus cut offsets. For example, `CCCGGG | 3|3` means SmaI cuts blunt after motif base 3 on both strands.",
                                );
                                ui.label("");
                                ui.end_row();
                                for suggestion in &displayed_suggestions {
                                    let suggestion_help =
                                        Self::gibson_destination_suggestion_help_text(suggestion);
                                    ui.add_sized(
                                        [SUGGESTION_COL_WIDTH, 0.0],
                                        egui::Label::new(
                                            egui::RichText::new(
                                                Self::gibson_destination_suggestion_display_label(
                                                    suggestion,
                                                ),
                                            )
                                            .small(),
                                        )
                                        .wrap(),
                                    )
                                    .on_hover_text(suggestion_help.clone());
                                    ui.add_sized(
                                        [ENDS_COL_WIDTH, 0.0],
                                        egui::Label::new(
                                            egui::RichText::new(
                                                Self::gibson_destination_suggestion_geometry_label(
                                                    suggestion,
                                                ),
                                            )
                                            .small(),
                                        )
                                        .wrap(),
                                    )
                                    .on_hover_text(suggestion_help.clone());
                                    ui.add_sized(
                                        [FEATURE_COL_WIDTH, 0.0],
                                        egui::Label::new(
                                            egui::RichText::new(
                                                Self::gibson_destination_suggestion_feature_label(
                                                    suggestion,
                                                ),
                                            )
                                            .small(),
                                        )
                                        .wrap(),
                                    )
                                    .on_hover_text(suggestion_help.clone());
                                    ui.add_sized(
                                        [CUT_COL_WIDTH, 0.0],
                                        egui::Label::new(
                                            egui::RichText::new(
                                                Self::gibson_destination_suggestion_cut_label(
                                                    suggestion,
                                                ),
                                            )
                                            .monospace()
                                            .small(),
                                        )
                                        .wrap(),
                                    )
                                    .on_hover_text(suggestion_help.clone());
                                    if ui
                                        .button("Use")
                                        .on_hover_text(format!(
                                            "{}\nFill the visible cut/opening edge fields from this suggestion.",
                                            suggestion_help
                                        ))
                                        .clicked()
                                    {
                                        selected_opening_suggestion = Some(suggestion.clone());
                                    }
                                    ui.end_row();
                                }
                            });
                    });
                if hidden_unique_cutters > 0 && !self.gibson_show_all_unique_cutters {
                    if ui
                        .button(format!("Search other unique cleavage sites ({hidden_unique_cutters})"))
                        .on_hover_text(
                            "Reveal unique cutters that are not named by the current MCS annotation.",
                        )
                        .clicked()
                    {
                        self.open_gibson_show_all_unique_cutters();
                    }
                } else if self.gibson_show_all_unique_cutters
                    && suggestions.iter().any(|row| row.in_mcs_context)
                {
                    if ui
                        .button("Show MCS-linked cutters only")
                        .on_hover_text(
                            "Collapse the list back to unique cutters named by the current MCS annotation.",
                        )
                        .clicked()
                    {
                        self.gibson_show_all_unique_cutters = false;
                    }
                }
            }
            Ok(_) => {
                if !self.gibson_destination_seq_id.trim().is_empty() {
                    ui.small(
                        "No unique-cutter suggestions were found for the current destination. You can still enter cut/opening edges directly or use the active selection.",
                    );
                    if ui
                        .button("Search unique cleavage sites")
                        .on_hover_text(
                            "Re-scan all known restriction enzymes for unique cleavage sites on the current destination.",
                        )
                        .clicked()
                    {
                        self.open_gibson_show_all_unique_cutters();
                        self.gibson_status = "No unique cleavage sites were found on the current destination after scanning all known restriction enzymes.".to_string();
                    }
                }
            }
            Err(err) => {
                ui.small(format!("Opening suggestions unavailable: {err}"));
            }
        }
        if let Some(suggestion) = selected_opening_suggestion.as_ref() {
            self.apply_gibson_destination_opening_suggestion(suggestion);
        }
        if let Some(mut opening_detail_text) = self.gibson_opening_detail_text(
            destination_opening_suggestions
                .as_ref()
                .ok()
                .map_or(&[], |rows| rows),
        ) {
            ui.small("Opening sketch");
            egui::Frame::group(ui.style()).show(ui, |ui| {
                ui.small("Exact destination sequence at the chosen cut/opening:");
                let opening_rows = opening_detail_text.lines().count().clamp(2, 4);
                ui.add(
                    egui::TextEdit::multiline(&mut opening_detail_text)
                        .desired_rows(opening_rows)
                        .desired_width(f32::INFINITY)
                        .font(egui::TextStyle::Monospace),
                );
                if let Some(mut opened_arm_text) = self.gibson_opened_destination_arm_text() {
                    ui.small("Destination after opening:");
                    let opened_rows = opened_arm_text.lines().count().clamp(2, 4);
                    ui.add(
                        egui::TextEdit::multiline(&mut opened_arm_text)
                            .desired_rows(opened_rows)
                            .desired_width(f32::INFINITY)
                            .font(egui::TextStyle::Monospace),
                    );
                }
                if let Some(preview) = self.gibson_preview_output.as_ref() {
                    if let Some(mut resolved_arm_text) =
                        self.gibson_resolved_destination_arm_text(preview)
                    {
                        ui.small("Resolved 5' destination ends after choosing the Gibson overlap:");
                        let resolved_rows = resolved_arm_text.lines().count().clamp(2, 4);
                        ui.add(
                            egui::TextEdit::multiline(&mut resolved_arm_text)
                                .desired_rows(resolved_rows)
                                .desired_width(f32::INFINITY)
                                .font(egui::TextStyle::Monospace),
                        );
                    }
                    if let Some(mut primer_construction_text) =
                        Self::gibson_primer_construction_text(preview)
                    {
                        ui.small("Insert primer construction:");
                        let primer_rows = primer_construction_text.lines().count().clamp(6, 10);
                        ui.add(
                            egui::TextEdit::multiline(&mut primer_construction_text)
                                .desired_rows(primer_rows)
                                .desired_width(f32::INFINITY)
                                .font(egui::TextStyle::Monospace),
                        );
                    }
                }
            });
        }

        ui.separator();
        ui.heading("Inserts");
        ui.small(
            "Insert rows are assembled in the listed order. One insert creates two terminal junctions; each additional insert adds one internal Gibson junction.",
        );
        ui.horizontal(|ui| {
            ui.label("1.");
            ui.label("sequence");
            egui::ComboBox::from_id_salt("gibson_insert_seq_id")
                .selected_text(if self.gibson_insert_seq_id.trim().is_empty() {
                    "(choose insert)"
                } else {
                    self.gibson_insert_seq_id.as_str()
                })
                .show_ui(ui, |ui| {
                    for seq_id in &sequence_ids {
                        if seq_id == &self.gibson_destination_seq_id {
                            continue;
                        }
                        ui.selectable_value(&mut self.gibson_insert_seq_id, seq_id.clone(), seq_id);
                    }
                });
            ui.text_edit_singleline(&mut self.gibson_insert_seq_id);
            ui.label("orientation");
            egui::ComboBox::from_id_salt("gibson_insert_orientation")
                .selected_text(self.gibson_insert_orientation.label())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.gibson_insert_orientation,
                        GibsonUiInsertOrientation::Forward,
                        GibsonUiInsertOrientation::Forward.label(),
                    );
                    ui.selectable_value(
                        &mut self.gibson_insert_orientation,
                        GibsonUiInsertOrientation::Reverse,
                        GibsonUiInsertOrientation::Reverse.label(),
                    );
                });
        });
        let mut remove_extra_index = None;
        let mut move_extra_up = None;
        let mut move_extra_down = None;
        let extra_insert_count = self.gibson_extra_inserts.len();
        for (idx, row) in self.gibson_extra_inserts.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                ui.label(format!("{}.", idx + 2));
                ui.label("sequence");
                egui::ComboBox::from_id_salt(format!("gibson_insert_seq_id_extra_{idx}"))
                    .selected_text(if row.seq_id.trim().is_empty() {
                        "(choose insert)"
                    } else {
                        row.seq_id.as_str()
                    })
                    .show_ui(ui, |ui| {
                        for seq_id in &sequence_ids {
                            if seq_id == &self.gibson_destination_seq_id {
                                continue;
                            }
                            ui.selectable_value(&mut row.seq_id, seq_id.clone(), seq_id);
                        }
                    });
                ui.text_edit_singleline(&mut row.seq_id);
                ui.label("orientation");
                egui::ComboBox::from_id_salt(format!("gibson_insert_orientation_extra_{idx}"))
                    .selected_text(row.orientation.label())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut row.orientation,
                            GibsonUiInsertOrientation::Forward,
                            GibsonUiInsertOrientation::Forward.label(),
                        );
                        ui.selectable_value(
                            &mut row.orientation,
                            GibsonUiInsertOrientation::Reverse,
                            GibsonUiInsertOrientation::Reverse.label(),
                        );
                    });
                if ui
                    .small_button("Up")
                    .on_hover_text("Move this insert earlier in the Gibson assembly order")
                    .clicked()
                {
                    move_extra_up = Some(idx);
                }
                if ui
                    .small_button("Down")
                    .on_hover_text("Move this insert later in the Gibson assembly order")
                    .clicked()
                    && idx + 1 < extra_insert_count
                {
                    move_extra_down = Some(idx);
                }
                if ui
                    .small_button("Remove")
                    .on_hover_text("Remove this insert row from the Gibson assembly order")
                    .clicked()
                {
                    remove_extra_index = Some(idx);
                }
            });
        }
        if let Some(idx) = move_extra_up {
            if idx == 0 {
                std::mem::swap(
                    &mut self.gibson_insert_seq_id,
                    &mut self.gibson_extra_inserts[0].seq_id,
                );
                std::mem::swap(
                    &mut self.gibson_insert_orientation,
                    &mut self.gibson_extra_inserts[0].orientation,
                );
            } else {
                self.gibson_extra_inserts.swap(idx, idx - 1);
            }
        }
        if let Some(idx) = move_extra_down {
            self.gibson_extra_inserts.swap(idx, idx + 1);
        }
        if let Some(idx) = remove_extra_index {
            self.gibson_extra_inserts.remove(idx);
        }
        if ui
            .button("Add Insert")
            .on_hover_text("Append another insert to the Gibson assembly order")
            .clicked()
        {
            self.gibson_extra_inserts.push(GibsonUiInsertRow::default());
        }
        self.refresh_gibson_output_id_hint_default();

        ui.separator();
        ui.heading("Design Targets");
        ui.horizontal(|ui| {
            ui.label("overlap min/max bp");
            ui.add(egui::TextEdit::singleline(&mut self.gibson_overlap_bp_min).desired_width(60.0));
            ui.add(egui::TextEdit::singleline(&mut self.gibson_overlap_bp_max).desired_width(60.0));
            Self::gibson_tm_label(
                ui,
                "minimum overlap ",
                " (°C)",
                egui::TextStyle::Body,
                false,
            )
            .on_hover_text(&tm_help);
            ui.add(
                egui::TextEdit::singleline(&mut self.gibson_minimum_overlap_tm_celsius)
                    .desired_width(80.0),
            );
        });
        ui.horizontal(|ui| {
            Self::gibson_tm_label(ui, "priming ", " window (°C)", egui::TextStyle::Body, false)
                .on_hover_text(&tm_help);
            ui.add(
                egui::TextEdit::singleline(&mut self.gibson_priming_tm_min_celsius)
                    .desired_width(80.0),
            );
            ui.add(
                egui::TextEdit::singleline(&mut self.gibson_priming_tm_max_celsius)
                    .desired_width(80.0),
            );
            ui.label("priming length bp");
            ui.add(
                egui::TextEdit::singleline(&mut self.gibson_priming_length_min_bp)
                    .desired_width(60.0),
            );
            ui.add(
                egui::TextEdit::singleline(&mut self.gibson_priming_length_max_bp)
                    .desired_width(60.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("new unique site").on_hover_text(
                "Optional REBASE enzyme name to introduce once on the assembled product through one engineered terminal Gibson overlap. Current v1 supports defined-site, single-insert plans and palindromic cutters.",
            );
            ui.add(
                egui::TextEdit::singleline(
                    &mut self.gibson_unique_restriction_site_enzyme_name,
                )
                .desired_width(120.0)
                .hint_text("e.g. EcoRI"),
            );
            if ui
                .small_button("Clear")
                .on_hover_text("Remove the optional unique-site engineering request")
                .clicked()
            {
                self.gibson_unique_restriction_site_enzyme_name.clear();
            }
        });
        ui.horizontal(|ui| {
            ui.label("output id hint");
            ui.add(
                egui::TextEdit::singleline(&mut self.gibson_output_id_hint)
                    .desired_width(f32::INFINITY),
            );
        });
        ui.small(tm_help.as_str());
        ui.small(
            "Optional unique-site request: ask Gibson to introduce one new unique restriction site, such as EcoRI, by mutating one terminal overlap if the resulting assembled product still stays uniquely cut there.",
        );

        ui.separator();
        Self::gibson_tm_label(ui, "", " Model", egui::TextStyle::Heading, true);
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.small("The Gibson specialist and CLI share one deterministic Tm model so preview numbers stay consistent across interfaces.");
            ui.small(tm_help.as_str());
        });

        ui.separator();
        ui.heading("Review");
        ui.horizontal(|ui| {
            if ui
                .add_enabled(preview_ready, egui::Button::new("Preview Gibson Plan"))
                .on_hover_text(
                    if preview_ready {
                        "Resolve overlaps, derive Gibson primer suggestions, validate, and refresh the cartoon preview"
                    } else {
                        preview_ready_reason
                            .as_deref()
                            .unwrap_or("Fill in the Gibson opening first.")
                    },
                )
                .clicked()
            {
                self.run_gibson_preview();
            }
            if ui
                .button("Clear Preview")
                .on_hover_text("Clear the current Gibson preview payload")
                .clicked()
            {
                self.gibson_preview_output = None;
                self.gibson_preview_svg_uri.clear();
                self.gibson_status = "Cleared Gibson preview output".to_string();
            }
        });
        if let Some(reason) = preview_ready_reason.as_deref() {
            ui.small(format!("Preview becomes available once: {reason}"));
        }

        let mut selected_design_adjustment: Option<GibsonSuggestedDesignAdjustment> = None;
        if let Some(preview) = self.gibson_preview_output.clone() {
            ui.small(format!(
                "preview schema: {} | can_execute={}",
                preview.schema, preview.can_execute
            ));
            if !preview.errors.is_empty() {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("blocking errors: {}", preview.errors.len()),
                );
            }
            if !preview.warnings.is_empty() {
                ui.colored_label(
                    egui::Color32::from_rgb(170, 120, 20),
                    format!("warnings: {}", preview.warnings.len()),
                );
            }
            if !preview.errors.is_empty()
                || !preview.warnings.is_empty()
                || !preview.notes.is_empty()
            {
                ui.small("Target review");
                egui::Frame::group(ui.style()).show(ui, |ui| {
                    for (color, row) in self.gibson_target_review_rows(&preview) {
                        ui.horizontal_wrapped(|ui| {
                            ui.colored_label(color, "-");
                            Self::gibson_rich_text_label(ui, &row, egui::TextStyle::Small, color);
                        });
                    }
                });
                if !preview.suggested_design_adjustments.is_empty() {
                    ui.small("Suggested next adjustments");
                    egui::Frame::group(ui.style()).show(ui, |ui| {
                        for adjustment in &preview.suggested_design_adjustments {
                            ui.horizontal_wrapped(|ui| {
                                if ui.button(&adjustment.label).clicked() {
                                    selected_design_adjustment = Some(adjustment.clone());
                                }
                                ui.small(&adjustment.summary);
                            });
                            ui.small(&adjustment.rationale);
                        }
                    });
                }
                if let Some(unique_site) = preview.unique_restriction_site.as_ref() {
                    ui.small("Requested unique site");
                    egui::Frame::group(ui.style()).show(ui, |ui| {
                        ui.small(&unique_site.message);
                        if unique_site.status.eq_ignore_ascii_case("engineered") {
                            ui.small(format!(
                                "{} terminal overlap '{}' now carries {} at overlap offset {} ({} mutated bp).",
                                unique_site.terminal_side,
                                unique_site.junction_id,
                                unique_site.recognition_sequence,
                                unique_site.motif_start_0based_in_overlap,
                                unique_site.mutated_bases
                            ));
                            ui.monospace(&unique_site.overlap_sequence);
                        }
                    });
                }
                ui.small("Preview findings");
                egui::Frame::group(ui.style()).show(ui, |ui| {
                    Self::render_gibson_preview_findings_rich(ui, &preview);
                });
                let mut findings_text = Self::gibson_preview_findings_text(&preview);
                let findings_rows = findings_text.lines().count().clamp(4, 10);
                ui.small("Preview findings (copyable)");
                ui.add(
                    egui::TextEdit::multiline(&mut findings_text)
                        .desired_rows(findings_rows)
                        .desired_width(f32::INFINITY)
                        .font(egui::TextStyle::Monospace),
                );
            }

            ui.separator();
            ui.label("Resolved junctions");
            egui::Grid::new("gibson_preview_junctions")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("junction");
                    ui.strong("members");
                    ui.strong("overlap");
                    Self::gibson_tm_label(ui, "", " (°C)", egui::TextStyle::Button, true)
                        .on_hover_text(&tm_help);
                    ui.end_row();
                    for junction in &preview.resolved_junctions {
                        ui.monospace(&junction.junction_id);
                        ui.label(format!(
                            "{} -> {}",
                            junction.left_member_id, junction.right_member_id
                        ));
                        ui.monospace(format!(
                            "{} bp | {}",
                            junction.overlap_bp, junction.overlap_sequence
                        ));
                        ui.label(format!("{:.1} °C", junction.overlap_tm_celsius))
                            .on_hover_text(&tm_help);
                        ui.end_row();
                    }
                });
            let mut junctions_text = Self::gibson_preview_junctions_text(&preview);
            if !junctions_text.trim().is_empty() {
                let junction_rows = junctions_text.lines().count().clamp(4, 10);
                ui.small("Resolved junctions (copyable)");
                ui.add(
                    egui::TextEdit::multiline(&mut junctions_text)
                        .desired_rows(junction_rows)
                        .desired_width(f32::INFINITY)
                        .font(egui::TextStyle::Monospace),
                );
            }

            ui.separator();
            ui.label("Primer suggestions");
            egui::Grid::new("gibson_preview_primers")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("primer");
                    ui.strong("full sequence");
                    ui.strong("5' overlap");
                    ui.strong("3' priming");
                    ui.end_row();
                    for primer in &preview.primer_suggestions {
                        ui.monospace(&primer.primer_id);
                        ui.monospace(&primer.full_sequence);
                        ui.small(format!(
                            "{} ({} bp, {:.1} °C)",
                            primer.overlap_5prime.sequence,
                            primer.overlap_5prime.length_bp,
                            primer.overlap_5prime.tm_celsius
                        ))
                        .on_hover_text(&tm_help);
                        ui.small(format!(
                            "{} ({} bp, {:.1} °C, hits={})",
                            primer.priming_3prime.sequence,
                            primer.priming_3prime.length_bp,
                            primer.priming_3prime.tm_celsius,
                            primer.priming_3prime.anneal_hits
                        ))
                        .on_hover_text(&tm_help);
                        ui.end_row();
                    }
                });
            let mut primers_text = Self::gibson_preview_primers_text(&preview);
            if !primers_text.trim().is_empty() {
                let primer_rows = primers_text.lines().count().clamp(4, 12);
                ui.small("Primer suggestions (copyable)");
                ui.add(
                    egui::TextEdit::multiline(&mut primers_text)
                        .desired_rows(primer_rows)
                        .desired_width(f32::INFINITY)
                        .font(egui::TextStyle::Monospace),
                );
            }

            ui.separator();
            ui.label("Cartoon preview");
            if !self.gibson_preview_svg_uri.trim().is_empty() {
                ui.add(
                    egui::Image::from_uri(self.gibson_preview_svg_uri.clone())
                        .max_width(ui.available_width())
                        .max_height(420.0)
                        .shrink_to_fit(),
                );
            } else {
                ui.small("No in-window cartoon preview is currently loaded.");
            }
        } else {
            ui.small("No Gibson preview yet. Use 'Preview Gibson Plan' to resolve junctions, primers, and cartoon bindings.");
        }
        if let Some(adjustment) = selected_design_adjustment.as_ref() {
            self.apply_gibson_design_adjustment(adjustment);
        }

        ui.separator();
        ui.heading("Outputs");
        if let Some(preview) = self.gibson_preview_output.as_ref() {
            let planned_product_id = if self.gibson_output_id_hint.trim().is_empty() {
                format!(
                    "{}_with_{}",
                    preview.destination.seq_id,
                    preview
                        .fragments
                        .iter()
                        .map(|fragment| fragment.seq_id.as_str())
                        .collect::<Vec<_>>()
                        .join("_")
                )
            } else {
                self.gibson_output_id_hint.trim().to_string()
            };
            let total_insert_bp = preview
                .fragments
                .iter()
                .map(|fragment| fragment.length_bp)
                .sum::<usize>();
            let planned_product_length_bp = if preview
                .destination
                .opening_mode
                .eq_ignore_ascii_case("existing_termini")
            {
                preview.destination.length_bp + total_insert_bp
            } else {
                preview
                    .destination
                    .length_bp
                    .saturating_sub(preview.destination.removed_span_bp.unwrap_or(0))
                    + total_insert_bp
            };
            let planned_product_topology = match self.build_gibson_plan_from_ui() {
                Ok(plan)
                    if plan
                        .product
                        .topology
                        .trim()
                        .eq_ignore_ascii_case("circular") =>
                {
                    "circular"
                }
                Ok(_) => "linear",
                Err(_) => preview.destination.actual_topology.as_str(),
            };
            ui.label("Planned output nodes");
            egui::Grid::new("gibson_planned_outputs")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("kind");
                    ui.strong("sequence id");
                    ui.strong("shape");
                    ui.end_row();
                    for primer in &preview.primer_suggestions {
                        let primer_kind = match primer.side.as_str() {
                            "left_insert_primer" => "Left insert primer",
                            "right_insert_primer" => "Right insert primer",
                            _ => primer.side.as_str(),
                        };
                        ui.label(primer_kind);
                        ui.monospace(&primer.primer_id);
                        ui.label(format!("{} bp, linear", primer.full_sequence.len()));
                        ui.end_row();
                    }
                    ui.label("Assembled product");
                    ui.monospace(planned_product_id);
                    ui.label(format!(
                        "{} bp, {}",
                        planned_product_length_bp, planned_product_topology
                    ));
                    ui.end_row();
                });
            ui.small(
                "Applying creates these primer/product sequence nodes through one Gibson operation. If an ID already exists, GENtle adds a numeric suffix deterministically.",
            );
        } else {
            ui.small(
                "Preview first to see the exact primer and assembled-product nodes that Gibson apply will create.",
            );
        }
        if self.gibson_is_multi_insert_draft()
            && self.gibson_opening_mode == GibsonUiOpeningMode::ExistingTermini
        {
            ui.colored_label(
                ui.visuals().warn_fg_color,
                Self::gibson_multi_insert_defined_opening_note(),
            );
        }
        ui.horizontal_wrapped(|ui| {
            let apply_enabled = self
                .gibson_preview_output
                .as_ref()
                .map(|preview| preview.can_execute)
                .unwrap_or(false);
            if ui
                .add_enabled(
                    apply_enabled,
                    egui::Button::new("Apply Gibson Cloning"),
                )
                .on_hover_text(
                    if apply_enabled {
                        "Create new sequence nodes for the Gibson primers and assembled product through one shared engine operation"
                    } else {
                        "Apply becomes available after Preview Gibson Plan reports no blocking errors."
                    },
                )
                .clicked()
            {
                self.run_gibson_apply();
            }
            if ui
                .add_enabled(preview_ready, egui::Button::new("Export Plan JSON..."))
                .on_hover_text("Export the canonical Gibson plan JSON currently implied by this window")
                .clicked()
            {
                self.export_gibson_plan_json();
            }
            if ui
                .add_enabled(
                    self.gibson_preview_output.is_some(),
                    egui::Button::new("Export Preview JSON..."),
                )
                .on_hover_text("Export the resolved Gibson preview payload")
                .clicked()
            {
                self.export_gibson_preview_json();
            }
            if ui
                .add_enabled(
                    self.gibson_preview_output.is_some(),
                    egui::Button::new("Export Primer Summary..."),
                )
                .on_hover_text("Export a text summary of the current Gibson primer suggestions")
                .clicked()
            {
                self.export_gibson_primer_summary();
            }
            if ui
                .add_enabled(
                    self.gibson_preview_output.is_some(),
                    egui::Button::new("Export Cartoon SVG..."),
                )
                .on_hover_text("Export the currently resolved Gibson cartoon as SVG")
                .clicked()
            {
                self.export_gibson_cartoon_svg();
            }
            if ui
                .add_enabled(
                    self.gibson_preview_output.is_some(),
                    egui::Button::new("Open in Routine Assistant"),
                )
                .on_hover_text("Hand the current Gibson preview into the existing routine workflow when possible")
                .clicked()
            {
                self.handoff_gibson_preview_to_routine_assistant();
            }
            if ui
                .button("Cancel")
                .on_hover_text("Close the Gibson specialist without applying anything")
                .clicked()
            {
                cancel_clicked = true;
            }
        });

        if !self.gibson_status.trim().is_empty() {
            ui.separator();
            ui.small("Status (copyable)");
            let mut status_text = Self::humanize_gibson_ui_text(&self.gibson_status);
            let status_rows = status_text.lines().count().clamp(2, 5);
            ui.add(
                egui::TextEdit::multiline(&mut status_text)
                    .desired_rows(status_rows)
                    .desired_width(f32::INFINITY)
                    .font(egui::TextStyle::Monospace),
            );
        }
        ui.separator();
        egui::CollapsingHeader::new("Quick Help")
            .default_open(true)
            .show(ui, |ui| {
                ui.small(
                    "GENtle starts with unique cutters named by the current MCS annotation because those are usually the most meaningful cloning openings.",
                );
                ui.small(
                    "Use a named cutter such as `SmaI` when you want one concrete restriction site and its end geometry.",
                );
                ui.small(
                    "The left/right cut edges are the primer-design anchors. For sticky cutters, inspect both 5' arm ends independently rather than treating them as one generic region.",
                );
                ui.small(
                    "Use `Show other unique cutters` when you want to look beyond the MCS and search the rest of the destination for single-cutters.",
                );
                ui.small(
                    "The Feature column shows the overlapping gene when available, otherwise another overlapping feature; MCS-linked cutters also show the MCS location.",
                );
                ui.small(
                    "The Cut column shows the REBASE recognition motif and the stored cut offsets used for that suggestion.",
                );
                ui.small(
                    "Equal left/right cut edges mean a blunt cutpoint. Different left/right edges mean the two primer-relevant termini of a sticky cut.",
                );
                ui.small(tm_help.as_str());
                ui.small(active_context_summary);
            });
        cancel_clicked
    }

    pub(super) fn render_gibson_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_gibson_dialog {
            return;
        }
        let mut open = self.show_gibson_dialog;
        let mut cancel_clicked = false;
        let viewport_id = Self::gibson_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Gibson",
            egui::Id::new(("hosted_gibson_window", viewport_id)),
            viewport_id,
            Vec2::new(1040.0, 820.0),
            Vec2::new(760.0, 560.0),
        );
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                    egui::ScrollArea::vertical()
                        .id_salt("gibson_embedded_scroll")
                        .auto_shrink([false, false])
                        .show(ui, |ui| cancel_clicked = self.render_gibson_contents(ui));
                });
            } else {
                crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                    egui::ScrollArea::vertical()
                        .id_salt("gibson_viewport_scroll")
                        .auto_shrink([false, false])
                        .show(ui, |ui| cancel_clicked = self.render_gibson_contents(ui));
                });
                if Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        if cancel_clicked {
            open = false;
        }
        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            open = false;
        }
        if !open {
            self.cancel_gibson_dialog();
        } else {
            self.show_gibson_dialog = true;
        }
    }
}
