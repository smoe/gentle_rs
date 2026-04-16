//! Variant-linked promoter follow-up expert for `MainAreaDna`.
//!
//! This keeps the GUI orchestration for promoter-SNP-to-luciferase follow-up
//! close to the sequence window while still routing all biology/business logic
//! through shared engine operations.

use super::*;

impl MainAreaDna {
    pub(super) fn feature_kind_supports_variant_followup(kind_upper: &str) -> bool {
        kind_upper == "VARIATION"
    }

    pub(super) fn feature_supports_variant_followup(&self, feature_id: usize) -> bool {
        self.dna
            .read()
            .ok()
            .and_then(|dna| {
                dna.features().get(feature_id).map(|feature| {
                    let kind_upper = feature.kind.to_string().trim().to_ascii_uppercase();
                    Self::feature_kind_supports_variant_followup(kind_upper.as_str())
                })
            })
            .unwrap_or(false)
    }

    pub(super) fn open_variant_followup_for_feature(
        &mut self,
        feature_id: usize,
        source: &str,
    ) -> bool {
        if !self.feature_supports_variant_followup(feature_id) {
            self.op_status = format!(
                "Could not open Variant Follow-up from {source}: selected feature is not a variation"
            );
            return false;
        }
        if let Err(err) = self.seed_variant_followup_defaults_for_feature(feature_id) {
            self.op_status = format!("Could not open Variant Follow-up from {source}: {err}");
            self.op_error_popup = Some(err);
            return false;
        }
        self.show_variant_followup_window = true;
        self.op_status = format!("Opened Variant Follow-up from {source}");
        self.op_error_popup = None;
        true
    }

    fn variant_followup_viewport_id(
        source_seq_id: &str,
        feature_id: Option<usize>,
    ) -> egui::ViewportId {
        egui::ViewportId::from_hash_of((
            "variant_followup_viewport",
            source_seq_id,
            feature_id.unwrap_or(usize::MAX),
        ))
    }

    fn variant_followup_window_title(ui: &VariantFollowupUiState) -> String {
        let label = ui.variant_label_or_id.trim();
        if label.is_empty() {
            format!("Variant Follow-up ({})", ui.source_seq_id)
        } else {
            format!("Variant Follow-up - {} ({})", label, ui.source_seq_id)
        }
    }

    fn variant_followup_help_text() -> &'static str {
        "ClawBio-facing promoter follow-up for a selected variation: derive transcript-based promoter windows, summarize promoter context, propose a luciferase reporter fragment, materialize matched alleles, and preview a mammalian promoterless luciferase reporter pair."
    }

    fn variant_followup_sequence_exists(&self, seq_id: &str) -> bool {
        if seq_id.trim().is_empty() {
            return false;
        }
        if self.seq_id.as_deref() == Some(seq_id) {
            return true;
        }
        self.engine
            .as_ref()
            .and_then(|engine| {
                engine
                    .read()
                    .ok()
                    .and_then(|guard| guard.state().sequences.get(seq_id).map(|_| ()))
            })
            .is_some()
    }

    fn variant_followup_feature_clone(
        &self,
        seq_id: &str,
        feature_id: usize,
    ) -> Option<gb_io::seq::Feature> {
        if self.seq_id.as_deref() == Some(seq_id) {
            return self.dna.read().ok()?.features().get(feature_id).cloned();
        }
        self.engine.as_ref().and_then(|engine| {
            engine.read().ok().and_then(|guard| {
                guard
                    .state()
                    .sequences
                    .get(seq_id)
                    .and_then(|dna| dna.features().get(feature_id).cloned())
            })
        })
    }

    fn variant_followup_extract_rsid_token(raw: &str) -> Option<String> {
        let mut current = String::new();
        let push_if_match = |token: &str| -> Option<String> {
            if token.len() > 2
                && token
                    .get(..2)
                    .map(|prefix| prefix.eq_ignore_ascii_case("rs"))
                    .unwrap_or(false)
                && token[2..].chars().all(|ch| ch.is_ascii_digit())
            {
                Some(format!("rs{}", &token[2..]))
            } else {
                None
            }
        };
        for ch in raw.chars() {
            if ch.is_ascii_alphanumeric() {
                current.push(ch);
            } else if let Some(found) = push_if_match(&current) {
                return Some(found);
            } else {
                current.clear();
            }
        }
        push_if_match(&current)
    }

    fn variant_followup_feature_variant_label(feature: &gb_io::seq::Feature) -> Option<String> {
        for key in [
            "db_xref",
            "label",
            "name",
            "variation",
            "standard_name",
            "note",
        ] {
            for value in feature.qualifier_values(key) {
                if let Some(rsid) = Self::variant_followup_extract_rsid_token(value) {
                    return Some(rsid);
                }
            }
        }
        Self::feature_tree_first_nonempty_qualifier(
            feature,
            &["label", "name", "variation", "standard_name", "db_xref"],
        )
    }

    fn variant_followup_suggested_token(raw: &str) -> String {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return "variant".to_string();
        }
        if let Some(rsid) = Self::variant_followup_extract_rsid_token(trimmed) {
            return Self::sanitize_export_name_component(&rsid, "variant");
        }
        Self::sanitize_export_name_component(trimmed, "variant")
    }

    fn seed_variant_followup_defaults_for_feature(
        &mut self,
        feature_id: usize,
    ) -> Result<(), String> {
        let source_seq_id = self
            .seq_id
            .clone()
            .ok_or_else(|| "No active sequence selected".to_string())?;
        let feature = self
            .variant_followup_feature_clone(&source_seq_id, feature_id)
            .ok_or_else(|| format!("Feature n-{feature_id} is no longer available"))?;
        let variant_label = Self::variant_followup_feature_variant_label(&feature)
            .unwrap_or_else(|| format!("variation_n{}", feature_id + 1));
        let gene_label =
            Self::feature_tree_first_nonempty_qualifier(&feature, &["gene"]).unwrap_or_default();
        let token = Self::variant_followup_suggested_token(&variant_label);
        self.variant_followup_ui = VariantFollowupUiState {
            source_seq_id,
            source_feature_id: Some(feature_id),
            variant_label_or_id: variant_label,
            gene_label,
            transcript_id: String::new(),
            promoter_upstream_bp: "1000".to_string(),
            promoter_downstream_bp: "200".to_string(),
            tfbs_focus_half_window_bp: "100".to_string(),
            retain_downstream_from_tss_bp: "200".to_string(),
            retain_upstream_beyond_variant_bp: "500".to_string(),
            max_candidates: "5".to_string(),
            fragment_output_id: format!("{token}_promoter_fragment"),
            reference_output_id: format!("{token}_promoter_reference"),
            alternate_output_id: format!("{token}_promoter_alternate"),
            reporter_backbone_seq_id: "gentle_mammalian_luciferase_backbone_v1".to_string(),
            reporter_backbone_path:
                "data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb".to_string(),
            reporter_output_prefix: format!("{token}_reporter"),
            cached_report: None,
            cached_candidates: None,
        };
        Ok(())
    }

    fn variant_followup_input_seq_id(&self) -> Result<String, String> {
        let seq_id = self.variant_followup_ui.source_seq_id.trim();
        if seq_id.is_empty() {
            Err("Variant Follow-up is not seeded from a source sequence yet".to_string())
        } else {
            Ok(seq_id.to_string())
        }
    }

    fn variant_followup_optional_text(value: &str) -> Option<String> {
        let trimmed = value.trim();
        (!trimmed.is_empty()).then(|| trimmed.to_string())
    }

    fn variant_followup_apply_report_defaults(&mut self, report: &VariantPromoterContextReport) {
        if self.variant_followup_ui.gene_label.trim().is_empty() {
            if let Some(gene_label) = report.chosen_gene_label.as_deref() {
                self.variant_followup_ui.gene_label = gene_label.to_string();
            }
        }
        if self.variant_followup_ui.transcript_id.trim().is_empty() {
            if let Some(transcript_id) = report.chosen_transcript_id.as_deref() {
                self.variant_followup_ui.transcript_id = transcript_id.to_string();
            }
        }
    }

    fn variant_followup_apply_candidate_defaults(
        &mut self,
        candidates: &PromoterReporterCandidateSet,
    ) {
        if self.variant_followup_ui.gene_label.trim().is_empty() {
            if let Some(gene_label) = candidates.chosen_gene_label.as_deref() {
                self.variant_followup_ui.gene_label = gene_label.to_string();
            }
        }
        if self.variant_followup_ui.transcript_id.trim().is_empty() {
            if let Some(transcript_id) = candidates.chosen_transcript_id.as_deref() {
                self.variant_followup_ui.transcript_id = transcript_id.to_string();
            }
        }
    }

    fn annotate_variant_followup_promoter_windows(&mut self) {
        let input = match self.variant_followup_input_seq_id() {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let upstream_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.promoter_upstream_bp,
            "promoter upstream bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let downstream_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.promoter_downstream_bp,
            "promoter downstream bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        self.variant_followup_ui.cached_report = None;
        self.variant_followup_ui.cached_candidates = None;
        self.apply_operation_with_feedback(Operation::AnnotatePromoterWindows {
            input,
            gene_label: Self::variant_followup_optional_text(&self.variant_followup_ui.gene_label),
            transcript_id: Self::variant_followup_optional_text(
                &self.variant_followup_ui.transcript_id,
            ),
            upstream_bp,
            downstream_bp,
            collapse_mode: PromoterWindowCollapseMode::Transcript,
        });
    }

    fn summarize_variant_followup_promoter_context(&mut self) {
        let input = match self.variant_followup_input_seq_id() {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let promoter_upstream_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.promoter_upstream_bp,
            "promoter upstream bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let promoter_downstream_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.promoter_downstream_bp,
            "promoter downstream bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let tfbs_focus_half_window_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.tfbs_focus_half_window_bp,
            "TFBS focus half-window bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let result = self.apply_operation_with_feedback_and_result(
            Operation::SummarizeVariantPromoterContext {
                input,
                variant_label_or_id: Self::variant_followup_optional_text(
                    &self.variant_followup_ui.variant_label_or_id,
                ),
                gene_label: Self::variant_followup_optional_text(
                    &self.variant_followup_ui.gene_label,
                ),
                transcript_id: Self::variant_followup_optional_text(
                    &self.variant_followup_ui.transcript_id,
                ),
                promoter_upstream_bp,
                promoter_downstream_bp,
                tfbs_focus_half_window_bp,
                path: None,
            },
        );
        if let Some(report) = result.and_then(|row| row.variant_promoter_context) {
            self.variant_followup_apply_report_defaults(&report);
            self.variant_followup_ui.cached_report = Some(report);
        }
    }

    fn suggest_variant_followup_reporter_fragments(&mut self) {
        let input = match self.variant_followup_input_seq_id() {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let retain_downstream_from_tss_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.retain_downstream_from_tss_bp,
            "retain downstream from TSS bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let retain_upstream_beyond_variant_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.retain_upstream_beyond_variant_bp,
            "retain upstream beyond variant bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let max_candidates = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.max_candidates,
            "max_candidates",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let result = self.apply_operation_with_feedback_and_result(
            Operation::SuggestPromoterReporterFragments {
                input,
                variant_label_or_id: Self::variant_followup_optional_text(
                    &self.variant_followup_ui.variant_label_or_id,
                ),
                gene_label: Self::variant_followup_optional_text(
                    &self.variant_followup_ui.gene_label,
                ),
                transcript_id: Self::variant_followup_optional_text(
                    &self.variant_followup_ui.transcript_id,
                ),
                retain_downstream_from_tss_bp,
                retain_upstream_beyond_variant_bp,
                max_candidates,
                path: None,
            },
        );
        if let Some(candidates) = result.and_then(|row| row.promoter_reporter_candidates) {
            self.variant_followup_apply_candidate_defaults(&candidates);
            self.variant_followup_ui.cached_candidates = Some(candidates);
        }
    }

    fn variant_followup_recommended_candidate(
        &self,
    ) -> Result<crate::engine::PromoterReporterFragmentCandidate, String> {
        let Some(candidates) = self.variant_followup_ui.cached_candidates.as_ref() else {
            return Err(
                "No reporter-fragment candidates are cached yet; run 'Propose reporter fragment' first"
                    .to_string(),
            );
        };
        candidates
            .candidates
            .iter()
            .find(|row| row.recommended)
            .cloned()
            .or_else(|| {
                candidates
                    .candidates
                    .iter()
                    .find(|row| row.candidate_id == candidates.recommended_candidate_id)
                    .cloned()
            })
            .ok_or_else(|| "No recommended promoter-reporter candidate is available".to_string())
    }

    fn extract_variant_followup_recommended_fragment(&mut self) {
        let input = match self.variant_followup_input_seq_id() {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let candidate = match self.variant_followup_recommended_candidate() {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let output_id =
            Self::variant_followup_optional_text(&self.variant_followup_ui.fragment_output_id);
        let result = self.apply_operation_with_feedback_and_result(Operation::ExtractRegion {
            input,
            from: candidate.start_0based,
            to: candidate.end_0based_exclusive,
            output_id,
        });
        if let Some(created_seq_id) = result
            .as_ref()
            .and_then(|row| row.created_seq_ids.first())
            .cloned()
        {
            self.variant_followup_ui.fragment_output_id = created_seq_id;
        }
    }

    fn materialize_variant_followup_alleles(&mut self) {
        let input = self
            .variant_followup_ui
            .fragment_output_id
            .trim()
            .to_string();
        if input.is_empty() {
            self.op_status =
                "Extract the recommended promoter fragment first so allele materialization has an input sequence"
                    .to_string();
            return;
        }
        if !self.variant_followup_sequence_exists(&input) {
            self.op_status = format!(
                "Promoter fragment '{}' is not available yet; extract the recommended fragment first",
                input
            );
            return;
        }
        let variant_label_or_id =
            Self::variant_followup_optional_text(&self.variant_followup_ui.variant_label_or_id);
        let reference_output_id =
            Self::variant_followup_optional_text(&self.variant_followup_ui.reference_output_id);
        let alternate_output_id =
            Self::variant_followup_optional_text(&self.variant_followup_ui.alternate_output_id);
        let reference_result =
            self.apply_operation_with_feedback_and_result(Operation::MaterializeVariantAllele {
                input: input.clone(),
                variant_label_or_id: variant_label_or_id.clone(),
                allele: VariantAlleleChoice::Reference,
                output_id: reference_output_id,
            });
        let Some(reference_seq_id) = reference_result
            .as_ref()
            .and_then(|row| row.created_seq_ids.first())
            .cloned()
        else {
            return;
        };
        self.variant_followup_ui.reference_output_id = reference_seq_id.clone();
        let alternate_result =
            self.apply_operation_with_feedback_and_result(Operation::MaterializeVariantAllele {
                input,
                variant_label_or_id,
                allele: VariantAlleleChoice::Alternate,
                output_id: alternate_output_id,
            });
        let Some(alternate_seq_id) = alternate_result
            .as_ref()
            .and_then(|row| row.created_seq_ids.first())
            .cloned()
        else {
            return;
        };
        self.variant_followup_ui.alternate_output_id = alternate_seq_id.clone();
        self.last_created_seq_ids = vec![reference_seq_id.clone(), alternate_seq_id.clone()];
        self.export_pool_inputs_text = self.last_created_seq_ids.join(", ");
        self.op_status = format!(
            "Materialized matched promoter alleles: '{}', '{}'",
            reference_seq_id, alternate_seq_id
        );
        self.op_error_popup = None;
    }

    fn ensure_variant_followup_reporter_backbone_loaded(&mut self) -> Result<(), String> {
        let seq_id = self
            .variant_followup_ui
            .reporter_backbone_seq_id
            .trim()
            .to_string();
        if seq_id.is_empty() {
            return Err("Reporter backbone sequence id is empty".to_string());
        }
        if self.variant_followup_sequence_exists(&seq_id) {
            return Ok(());
        }
        let path = self
            .variant_followup_ui
            .reporter_backbone_path
            .trim()
            .to_string();
        if path.is_empty() {
            return Err("Reporter backbone path is empty".to_string());
        }
        let current_seq_id = self.seq_id.clone();
        let result = self.apply_operation_with_feedback_and_result(Operation::LoadFile {
            path,
            as_id: Some(seq_id.clone()),
        });
        if result.is_none() {
            return Err(self.op_status.clone());
        }
        if let Some(previous_seq_id) = current_seq_id.filter(|value| value != &seq_id) {
            let status = self.op_status.clone();
            let popup = self.op_error_popup.clone();
            let _ = self.open_sequence_by_id(&previous_seq_id);
            self.op_status = status;
            self.op_error_popup = popup;
        }
        Ok(())
    }

    fn preview_variant_followup_reporter_pair(&mut self) {
        let reference_fragment_seq_id = self
            .variant_followup_ui
            .reference_output_id
            .trim()
            .to_string();
        let alternate_fragment_seq_id = self
            .variant_followup_ui
            .alternate_output_id
            .trim()
            .to_string();
        if reference_fragment_seq_id.is_empty()
            || alternate_fragment_seq_id.is_empty()
            || !self.variant_followup_sequence_exists(&reference_fragment_seq_id)
            || !self.variant_followup_sequence_exists(&alternate_fragment_seq_id)
        {
            self.materialize_variant_followup_alleles();
        }
        let reference_fragment_seq_id = self
            .variant_followup_ui
            .reference_output_id
            .trim()
            .to_string();
        let alternate_fragment_seq_id = self
            .variant_followup_ui
            .alternate_output_id
            .trim()
            .to_string();
        if reference_fragment_seq_id.is_empty()
            || alternate_fragment_seq_id.is_empty()
            || !self.variant_followup_sequence_exists(&reference_fragment_seq_id)
            || !self.variant_followup_sequence_exists(&alternate_fragment_seq_id)
        {
            return;
        }
        if let Err(err) = self.ensure_variant_followup_reporter_backbone_loaded() {
            self.op_status = err.clone();
            self.op_error_popup = Some(err);
            return;
        }
        let backbone_seq_id = self
            .variant_followup_ui
            .reporter_backbone_seq_id
            .trim()
            .to_string();
        let prefix = {
            let preferred = self.variant_followup_ui.reporter_output_prefix.trim();
            if preferred.is_empty() {
                Self::variant_followup_suggested_token(
                    &self.variant_followup_ui.variant_label_or_id,
                )
            } else {
                Self::sanitize_export_name_component(preferred, "reporter_pair")
            }
        };
        self.variant_followup_ui.reporter_output_prefix = prefix.clone();
        let reference_assembly_prefix = format!("{prefix}_reference_assembly");
        let alternate_assembly_prefix = format!("{prefix}_alternate_assembly");
        let reference_reporter_id = format!("{prefix}_reference");
        let alternate_reporter_id = format!("{prefix}_alternate");

        let reference_ligation =
            self.apply_operation_with_feedback_and_result(Operation::Ligation {
                inputs: vec![reference_fragment_seq_id, backbone_seq_id.clone()],
                circularize_if_possible: false,
                output_id: None,
                protocol: LigationProtocol::Blunt,
                output_prefix: Some(reference_assembly_prefix),
                unique: Some(false),
            });
        let Some(reference_assembly_id) = reference_ligation
            .as_ref()
            .and_then(|row| row.created_seq_ids.first())
            .cloned()
        else {
            return;
        };
        let reference_branch = self.apply_operation_with_feedback_and_result(Operation::Branch {
            input: reference_assembly_id,
            output_id: Some(reference_reporter_id.clone()),
        });
        let Some(reference_preview_id) = reference_branch
            .as_ref()
            .and_then(|row| row.created_seq_ids.first())
            .cloned()
        else {
            return;
        };

        let alternate_ligation =
            self.apply_operation_with_feedback_and_result(Operation::Ligation {
                inputs: vec![alternate_fragment_seq_id, backbone_seq_id],
                circularize_if_possible: false,
                output_id: None,
                protocol: LigationProtocol::Blunt,
                output_prefix: Some(alternate_assembly_prefix),
                unique: Some(false),
            });
        let Some(alternate_assembly_id) = alternate_ligation
            .as_ref()
            .and_then(|row| row.created_seq_ids.first())
            .cloned()
        else {
            return;
        };
        let alternate_branch = self.apply_operation_with_feedback_and_result(Operation::Branch {
            input: alternate_assembly_id,
            output_id: Some(alternate_reporter_id.clone()),
        });
        let Some(alternate_preview_id) = alternate_branch
            .as_ref()
            .and_then(|row| row.created_seq_ids.first())
            .cloned()
        else {
            return;
        };
        self.last_created_seq_ids =
            vec![reference_preview_id.clone(), alternate_preview_id.clone()];
        self.export_pool_inputs_text = self.last_created_seq_ids.join(", ");
        self.op_status = format!(
            "Previewed luciferase reporter pair: '{}', '{}'",
            reference_preview_id, alternate_preview_id
        );
        self.op_error_popup = None;
    }

    fn render_variant_followup_report_summary(&mut self, ui: &mut egui::Ui) {
        let Some(report) = self.variant_followup_ui.cached_report.as_ref() else {
            ui.small(
                egui::RichText::new(
                    "No promoter-context summary cached yet. Run 'Summarize promoter context' to capture the ClawBio handoff record.",
                )
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
            return;
        };
        ui.group(|ui| {
            ui.label(egui::RichText::new("Promoter context summary").strong());
            egui::Grid::new((
                "variant_followup_report_grid",
                report.seq_id.as_str(),
                report.variant_label.as_str(),
            ))
            .num_columns(2)
            .spacing([12.0, 4.0])
            .show(ui, |ui| {
                ui.small("Variant");
                ui.monospace(&report.variant_label);
                ui.end_row();
                ui.small("Chosen gene");
                ui.label(
                    report
                        .chosen_gene_label
                        .clone()
                        .unwrap_or_else(|| "<auto>".to_string()),
                );
                ui.end_row();
                ui.small("Chosen transcript");
                ui.label(
                    report
                        .chosen_transcript_label
                        .clone()
                        .or_else(|| report.chosen_transcript_id.clone())
                        .unwrap_or_else(|| "<auto>".to_string()),
                );
                ui.end_row();
                ui.small("Promoter overlap");
                ui.label(if report.promoter_overlap { "yes" } else { "no" });
                ui.end_row();
                ui.small("Signed TSS distance");
                ui.label(
                    report
                        .signed_tss_distance_bp
                        .map(|value| format!("{value} bp"))
                        .unwrap_or_else(|| "n/a".to_string()),
                );
                ui.end_row();
                ui.small("Transcript status");
                ui.label(&report.transcript_ambiguity_status);
                ui.end_row();
                ui.small("TFBS near variant");
                ui.label(&report.tfbs_near_variant_status);
                ui.end_row();
                ui.small("Suggested assay");
                ui.label(if report.suggested_assay_ids.is_empty() {
                    "none".to_string()
                } else {
                    report.suggested_assay_ids.join(", ")
                });
                ui.end_row();
            });
            ui.small(
                egui::RichText::new(&report.rationale).color(egui::Color32::from_rgb(71, 85, 105)),
            );
            if let Some(tfbs_summary) = report.tfbs_region_summary.as_ref() {
                if !tfbs_summary.rows.is_empty() {
                    ui.add_space(6.0);
                    ui.small(
                        egui::RichText::new("TFBS near-variant factors")
                            .strong()
                            .color(egui::Color32::from_rgb(51, 65, 85)),
                    );
                    egui::Grid::new((
                        "variant_followup_tfbs_grid",
                        report.seq_id.as_str(),
                        report.variant_label.as_str(),
                    ))
                    .num_columns(4)
                    .spacing([10.0, 2.0])
                    .show(ui, |ui| {
                        ui.small("TF");
                        ui.small("focus");
                        ui.small("context");
                        ui.small("share");
                        ui.end_row();
                        for row in tfbs_summary.rows.iter().take(6) {
                            ui.label(&row.tf_name);
                            ui.monospace(row.focus_occurrences.to_string());
                            ui.monospace(row.context_occurrences.to_string());
                            ui.monospace(format!("{:.2}", row.focus_share_of_context_occurrences));
                            ui.end_row();
                        }
                    });
                }
            } else {
                ui.small(
                    egui::RichText::new(
                        "No TFBS annotations were available on this sequence yet; the summary records that explicitly rather than assuming no TFBS context exists.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
            }
        });
    }

    fn render_variant_followup_candidate_summary(&mut self, ui: &mut egui::Ui) {
        let Some(candidates) = self.variant_followup_ui.cached_candidates.as_ref() else {
            ui.small(
                egui::RichText::new(
                    "No reporter-fragment candidates cached yet. Run 'Propose reporter fragment' to rank transcript-aware promoter inserts.",
                )
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
            return;
        };
        ui.group(|ui| {
            ui.label(egui::RichText::new("Reporter-fragment candidates").strong());
            ui.small(
                egui::RichText::new(format!(
                    "Recommended candidate: {}",
                    candidates.recommended_candidate_id
                ))
                .color(egui::Color32::from_rgb(71, 85, 105)),
            );
            egui::Grid::new((
                "variant_followup_candidate_grid",
                candidates.seq_id.as_str(),
                candidates.recommended_candidate_id.as_str(),
            ))
            .num_columns(7)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.small("Rank");
                ui.small("Tx");
                ui.small("Strand");
                ui.small("Interval");
                ui.small("Length");
                ui.small("Rec");
                ui.small("Why");
                ui.end_row();
                for row in candidates.candidates.iter().take(5) {
                    ui.monospace(row.rank.to_string());
                    ui.label(&row.transcript_label);
                    ui.monospace(&row.strand);
                    ui.monospace(format!(
                        "{}..{}",
                        row.start_0based, row.end_0based_exclusive
                    ));
                    ui.monospace(format!("{} bp", row.length_bp));
                    ui.label(if row.recommended { "yes" } else { "" });
                    ui.small(&row.rationale);
                    ui.end_row();
                }
            });
        });
    }

    fn render_variant_followup_window_contents(&mut self, ui: &mut egui::Ui) {
        let engine_available = self.engine.is_some();
        let source_seq_id = self.variant_followup_ui.source_seq_id.clone();
        let current_seq_id = self.seq_id.clone().unwrap_or_default();
        let source_missing = !self.variant_followup_sequence_exists(&source_seq_id);
        let mut summary_params_changed = false;
        let mut candidate_params_changed = false;

        ui.horizontal_wrapped(|ui| {
            ui.label(
                egui::RichText::new("Window guide [?]")
                    .size(9.0)
                    .color(egui::Color32::from_rgb(71, 85, 105)),
            )
            .on_hover_text(Self::variant_followup_help_text());
            ui.label(
                egui::RichText::new(
                    "Generated promoter windows are transcript-derived evidence owned by the engine; this window only orchestrates the shared operations and displays the portable records.",
                )
                .size(9.0)
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
        });

        if source_missing {
            ui.colored_label(
                egui::Color32::from_rgb(180, 83, 9),
                format!(
                    "Source context sequence '{}' is not currently available in the engine state.",
                    source_seq_id
                ),
            );
        } else {
            ui.small(
                egui::RichText::new(format!(
                    "Source context: {}{}",
                    source_seq_id,
                    if current_seq_id.is_empty() || current_seq_id == source_seq_id {
                        String::new()
                    } else {
                        format!(" | current window: {current_seq_id}")
                    }
                ))
                .color(egui::Color32::from_rgb(71, 85, 105)),
            );
        }

        ui.separator();
        egui::Grid::new((
            "variant_followup_controls_grid",
            self.variant_followup_ui.source_seq_id.as_str(),
        ))
        .num_columns(2)
        .spacing([12.0, 6.0])
        .show(ui, |ui| {
            ui.label("Variant label / id");
            if ui
                .text_edit_singleline(&mut self.variant_followup_ui.variant_label_or_id)
                .changed()
            {
                summary_params_changed = true;
                candidate_params_changed = true;
            }
            ui.end_row();

            ui.label("Gene label");
            if ui
                .text_edit_singleline(&mut self.variant_followup_ui.gene_label)
                .changed()
            {
                summary_params_changed = true;
                candidate_params_changed = true;
            }
            ui.end_row();

            ui.label("Transcript id");
            if ui
                .text_edit_singleline(&mut self.variant_followup_ui.transcript_id)
                .changed()
            {
                summary_params_changed = true;
                candidate_params_changed = true;
            }
            ui.end_row();

            ui.label("Promoter window (up/down bp)");
            ui.horizontal_wrapped(|ui| {
                if ui
                    .text_edit_singleline(&mut self.variant_followup_ui.promoter_upstream_bp)
                    .changed()
                {
                    summary_params_changed = true;
                }
                ui.label("/");
                if ui
                    .text_edit_singleline(&mut self.variant_followup_ui.promoter_downstream_bp)
                    .changed()
                {
                    summary_params_changed = true;
                }
            });
            ui.end_row();

            ui.label("TFBS focus half-window bp");
            if ui
                .text_edit_singleline(&mut self.variant_followup_ui.tfbs_focus_half_window_bp)
                .changed()
            {
                summary_params_changed = true;
            }
            ui.end_row();

            ui.label("Reporter fragment (+downstream / +beyond SNP)");
            ui.horizontal_wrapped(|ui| {
                if ui
                    .text_edit_singleline(
                        &mut self.variant_followup_ui.retain_downstream_from_tss_bp,
                    )
                    .changed()
                {
                    candidate_params_changed = true;
                }
                ui.label("/");
                if ui
                    .text_edit_singleline(
                        &mut self.variant_followup_ui.retain_upstream_beyond_variant_bp,
                    )
                    .changed()
                {
                    candidate_params_changed = true;
                }
            });
            ui.end_row();

            ui.label("Max candidates");
            if ui
                .text_edit_singleline(&mut self.variant_followup_ui.max_candidates)
                .changed()
            {
                candidate_params_changed = true;
            }
            ui.end_row();

            ui.label("Fragment output id");
            ui.text_edit_singleline(&mut self.variant_followup_ui.fragment_output_id);
            ui.end_row();

            ui.label("Reference / alternate output ids");
            ui.horizontal_wrapped(|ui| {
                ui.text_edit_singleline(&mut self.variant_followup_ui.reference_output_id);
                ui.label("/");
                ui.text_edit_singleline(&mut self.variant_followup_ui.alternate_output_id);
            });
            ui.end_row();

            ui.label("Reporter backbone");
            ui.horizontal_wrapped(|ui| {
                ui.text_edit_singleline(&mut self.variant_followup_ui.reporter_backbone_seq_id);
                ui.label("@");
                ui.text_edit_singleline(&mut self.variant_followup_ui.reporter_backbone_path);
            });
            ui.end_row();

            ui.label("Reporter output prefix");
            ui.text_edit_singleline(&mut self.variant_followup_ui.reporter_output_prefix);
            ui.end_row();
        });

        if summary_params_changed {
            self.variant_followup_ui.cached_report = None;
        }
        if summary_params_changed || candidate_params_changed {
            self.variant_followup_ui.cached_candidates = None;
        }

        ui.separator();
        ui.horizontal_wrapped(|ui| {
            if ui
                .add_enabled(
                    engine_available && !source_missing,
                    egui::Button::new("Annotate promoter windows"),
                )
                .on_hover_text(
                    "Generate transcript-derived promoter features on the source sequence so promoter-proximal SNP reasoning does not depend on imported promoter annotations.",
                )
                .clicked()
            {
                self.annotate_variant_followup_promoter_windows();
            }
            if ui
                .add_enabled(
                    engine_available && !source_missing,
                    egui::Button::new("Summarize promoter context"),
                )
                .on_hover_text(
                    "Build the portable variant-promoter context record for ClawBio/OpenClaw handoff.",
                )
                .clicked()
            {
                self.summarize_variant_followup_promoter_context();
            }
            if ui
                .add_enabled(
                    engine_available && !source_missing,
                    egui::Button::new("Propose reporter fragment"),
                )
                .on_hover_text(
                    "Rank transcript-aware promoter fragments that retain TSS-proximal context and keep the variant inside the insert.",
                )
                .clicked()
            {
                self.suggest_variant_followup_reporter_fragments();
            }
        });

        ui.horizontal_wrapped(|ui| {
            let has_candidates = self.variant_followup_ui.cached_candidates.is_some();
            if ui
                .add_enabled(has_candidates, egui::Button::new("Extract recommended fragment"))
                .on_hover_text(
                    "Materialize the recommended promoter fragment through the shared ExtractRegion operation.",
                )
                .clicked()
            {
                self.extract_variant_followup_recommended_fragment();
            }
            let has_fragment = !self.variant_followup_ui.fragment_output_id.trim().is_empty();
            if ui
                .add_enabled(has_fragment, egui::Button::new("Make reference/alternate inserts"))
                .on_hover_text(
                    "Build matched SNV-specific promoter inserts from the extracted fragment.",
                )
                .clicked()
            {
                self.materialize_variant_followup_alleles();
            }
            let has_materialized = !self.variant_followup_ui.reference_output_id.trim().is_empty()
                && !self.variant_followup_ui.alternate_output_id.trim().is_empty();
            if ui
                .add_enabled(has_materialized, egui::Button::new("Preview luciferase pair"))
                .on_hover_text(
                    "Load the pinned local mammalian promoterless luciferase backbone if needed, then build reference/alternate reporter previews.",
                )
                .clicked()
            {
                self.preview_variant_followup_reporter_pair();
            }
        });

        ui.separator();
        self.render_variant_followup_report_summary(ui);
        ui.add_space(8.0);
        self.render_variant_followup_candidate_summary(ui);
    }

    pub(super) fn render_variant_followup_window(&mut self, ctx: &egui::Context) {
        if !self.show_variant_followup_window {
            return;
        }
        if self.variant_followup_ui.source_seq_id.trim().is_empty() {
            self.show_variant_followup_window = false;
            return;
        }
        let title = Self::variant_followup_window_title(&self.variant_followup_ui);
        let viewport_id = Self::variant_followup_viewport_id(
            &self.variant_followup_ui.source_seq_id,
            self.variant_followup_ui.source_feature_id,
        );
        let default_size = Vec2::new(980.0, 720.0);
        let min_size = Vec2::new(760.0, 520.0);
        let content_min_size = Vec2::new(900.0, 640.0);
        let builder = egui::ViewportBuilder::default()
            .with_title(title.clone())
            .with_inner_size([default_size.x, default_size.y])
            .with_min_inner_size([min_size.x, min_size.y]);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut open = self.show_variant_followup_window;
                egui::Window::new(title.clone())
                    .id(egui::Id::new(format!(
                        "variant_followup_window_embedded_{}_{}",
                        self.variant_followup_ui.source_seq_id,
                        self.variant_followup_ui
                            .source_feature_id
                            .map(|value| value.to_string())
                            .unwrap_or_else(|| "none".to_string())
                    )))
                    .open(&mut open)
                    .resizable(true)
                    .default_size(default_size)
                    .show(ctx, |ui| {
                        let backdrop_settings = current_window_backdrop_settings();
                        paint_window_backdrop(ui, WindowBackdropKind::Splicing, &backdrop_settings);
                        egui::ScrollArea::both()
                            .id_salt(format!(
                                "variant_followup_scroll_embedded_{}_{}",
                                self.variant_followup_ui.source_seq_id,
                                self.variant_followup_ui
                                    .source_feature_id
                                    .map(|value| value.to_string())
                                    .unwrap_or_else(|| "none".to_string())
                            ))
                            .auto_shrink([false, false])
                            .show(ui, |ui| {
                                scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                    ui,
                                    scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                                );
                                ui.set_min_size(content_min_size);
                                self.render_variant_followup_window_contents(ui);
                            });
                    });
                self.show_variant_followup_window = open;
                return;
            }

            crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                let backdrop_settings = current_window_backdrop_settings();
                paint_window_backdrop(ui, WindowBackdropKind::Splicing, &backdrop_settings);
                egui::ScrollArea::both()
                    .id_salt(format!(
                        "variant_followup_scroll_viewport_{}_{}",
                        self.variant_followup_ui.source_seq_id,
                        self.variant_followup_ui
                            .source_feature_id
                            .map(|value| value.to_string())
                            .unwrap_or_else(|| "none".to_string())
                    ))
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        scroll_input_policy::apply_scrollarea_keyboard_navigation(
                            ui,
                            scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                        );
                        ui.set_min_size(content_min_size);
                        self.render_variant_followup_window_contents(ui);
                    });
            });

            if crate::app::GENtleApp::viewport_close_requested_or_shortcut(ctx) {
                self.show_variant_followup_window = false;
            }
        });
    }
}
