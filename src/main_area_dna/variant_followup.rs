//! Variant-linked promoter follow-up expert for `MainAreaDna`.
//!
//! This keeps the GUI orchestration for promoter-SNP-to-luciferase follow-up
//! close to the sequence window while still routing all biology/business logic
//! through shared engine operations.

use super::*;
use crate::engine::{TfbsScoreTrackCorrelationMetric, TfbsScoreTrackCorrelationSignalSource};

#[derive(Clone, Debug)]
pub(super) struct VariantFollowupBundleArtifacts {
    pub(super) bundle_id: String,
    pub(super) variant_promoter_context_json: String,
    pub(super) promoter_reporter_candidates_json: String,
    pub(super) promoter_context_svg: String,
    pub(super) reference_reporter_svg: String,
    pub(super) alternate_reporter_svg: String,
    pub(super) report_md: String,
    pub(super) result_json: String,
    pub(super) commands_sh: String,
}

impl MainAreaDna {
    pub(super) fn feature_kind_supports_variant_followup(kind_upper: &str) -> bool {
        matches!(kind_upper, "VARIATION" | "GENE" | "MRNA" | "PROMOTER")
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
                "Could not open Promoter design from {source}: selected feature is not promoter-design relevant"
            );
            return false;
        }
        if let Err(err) = self.seed_variant_followup_defaults_for_feature(feature_id) {
            self.op_status = format!("Could not open Promoter design from {source}: {err}");
            self.op_error_popup = Some(err);
            return false;
        }
        self.log_promoter_design_status(
            "open requested",
            self.variant_followup_window_pending_initial_render,
        );
        self.variant_followup_window_pending_initial_render = true;
        self.show_variant_followup_window = true;
        self.log_promoter_design_status("window state stored", true);
        self.op_status = format!("Opened Promoter design from {source}");
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
        let label = ui.gene_label.trim();
        if !label.is_empty() {
            format!("Promoter design - {} ({})", label, ui.source_seq_id)
        } else {
            let label = ui.variant_label_or_id.trim();
            if label.is_empty() {
                format!("Promoter design ({})", ui.source_seq_id)
            } else {
                format!("Promoter design - {} ({})", label, ui.source_seq_id)
            }
        }
    }

    fn variant_followup_help_text() -> &'static str {
        "Promoter-focused design workspace for a selected sequence feature: derive transcript-based promoter windows, inspect positive-only TF motif score tracks, summarize promoter context, propose luciferase reporter fragments, materialize matched alleles when a variant is present, and preview a mammalian promoterless luciferase reporter pair."
    }

    fn promoter_design_default_motifs(gene_label: &str) -> String {
        match gene_label.trim().to_ascii_uppercase().as_str() {
            "TERT" => "SP1".to_string(),
            "TP73" => "TP73,SP1,BACH2,PATZ1".to_string(),
            _ => "SP1".to_string(),
        }
    }

    fn variant_followup_sequence_len(&self, seq_id: &str) -> Option<usize> {
        if seq_id.trim().is_empty() {
            return None;
        }
        if self.seq_id.as_deref() == Some(seq_id) {
            return self.dna.read().ok().map(|dna| dna.len());
        }
        self.engine.as_ref().and_then(|engine| {
            engine.read().ok().and_then(|guard| {
                guard
                    .state()
                    .sequences
                    .get(seq_id)
                    .map(crate::dna_sequence::DNAsequence::len)
            })
        })
    }

    fn variant_followup_has_variant_seed(&self) -> bool {
        !self
            .variant_followup_ui
            .variant_label_or_id
            .trim()
            .is_empty()
    }

    fn feature_seed_gene_label(feature: &gb_io::seq::Feature) -> String {
        Self::feature_tree_first_nonempty_qualifier(
            feature,
            &[
                "gene",
                "gene_name",
                "label",
                "name",
                "standard_name",
                "locus_tag",
            ],
        )
        .unwrap_or_default()
    }

    fn feature_seed_transcript_id(feature: &gb_io::seq::Feature) -> String {
        Self::feature_tree_first_nonempty_qualifier(
            feature,
            &["transcript_id", "label", "name", "standard_name"],
        )
        .unwrap_or_default()
    }

    fn feature_seed_token(
        feature: &gb_io::seq::Feature,
        variant_label: &str,
        gene_label: &str,
    ) -> String {
        if !variant_label.trim().is_empty() {
            return Self::variant_followup_suggested_token(variant_label);
        }
        if !gene_label.trim().is_empty() {
            return Self::sanitize_export_name_component(gene_label, "promoter");
        }
        Self::sanitize_export_name_component(&feature.kind.to_string(), "promoter")
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
        let feature_kind = feature.kind.to_string().trim().to_ascii_uppercase();
        let variant_label = if feature_kind == "VARIATION" {
            Self::variant_followup_feature_variant_label(&feature)
                .unwrap_or_else(|| format!("variation_n{}", feature_id + 1))
        } else {
            String::new()
        };
        let gene_label = Self::feature_seed_gene_label(&feature);
        let transcript_id = if feature_kind == "MRNA" {
            Self::feature_seed_transcript_id(&feature)
        } else {
            String::new()
        };
        let token = Self::feature_seed_token(&feature, &variant_label, &gene_label);
        let seq_len = self
            .variant_followup_sequence_len(&source_seq_id)
            .unwrap_or_default();
        self.variant_followup_ui = VariantFollowupUiState {
            source_seq_id,
            source_feature_id: Some(feature_id),
            variant_label_or_id: variant_label,
            gene_label: gene_label.clone(),
            transcript_id,
            score_track_motifs: Self::promoter_design_default_motifs(&gene_label),
            score_track_start_0based: "0".to_string(),
            score_track_end_0based_exclusive: seq_len.to_string(),
            score_track_value_kind: TfbsScoreTrackValueKind::LlrBits,
            score_track_clip_negative: true,
            score_track_correlation_metric: TfbsScoreTrackCorrelationMetric::Pearson,
            score_track_correlation_signal_source:
                TfbsScoreTrackCorrelationSignalSource::MaxStrands,
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
            cached_score_tracks: None,
            cached_report: None,
            cached_candidates: None,
        };
        Ok(())
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

    fn variant_followup_input_seq_id(&self) -> Result<String, String> {
        let seq_id = self.variant_followup_ui.source_seq_id.trim();
        if seq_id.is_empty() {
            Err("Promoter design is not seeded from a source sequence yet".to_string())
        } else {
            Ok(seq_id.to_string())
        }
    }

    fn variant_followup_optional_text(value: &str) -> Option<String> {
        let trimmed = value.trim();
        (!trimmed.is_empty()).then(|| trimmed.to_string())
    }

    fn promoter_design_parse_motif_tokens(raw: &str) -> Vec<String> {
        let mut seen = HashSet::new();
        let mut out = vec![];
        for token in raw
            .split(|ch: char| ch == ',' || ch == ';' || ch.is_whitespace())
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            let normalized = token.to_ascii_uppercase();
            if seen.insert(normalized) {
                out.push(token.to_string());
            }
        }
        out
    }

    fn variant_followup_bundle_stem(&self) -> String {
        let gene = Self::variant_followup_optional_text(&self.variant_followup_ui.gene_label)
            .map(|value| Self::sanitize_export_name_component(&value, "gene"))
            .unwrap_or_else(|| "gene".to_string());
        let variant = Self::sanitize_export_name_component(
            &self.variant_followup_ui.variant_label_or_id,
            "variant",
        );
        format!("{gene}_{variant}")
    }

    pub(super) fn variant_followup_bundle_artifacts(
        &self,
        reference_reporter_seq_id: &str,
        alternate_reporter_seq_id: &str,
    ) -> VariantFollowupBundleArtifacts {
        let stem = self.variant_followup_bundle_stem();
        VariantFollowupBundleArtifacts {
            bundle_id: format!("{stem}_promoter_reporter_handoff"),
            variant_promoter_context_json: "variant_promoter_context.json".to_string(),
            promoter_reporter_candidates_json: "promoter_reporter_candidates.json".to_string(),
            promoter_context_svg: format!("{stem}_promoter_context.svg"),
            reference_reporter_svg: format!(
                "{}.svg",
                Self::sanitize_export_name_component(
                    reference_reporter_seq_id,
                    "reporter_reference"
                )
            ),
            alternate_reporter_svg: format!(
                "{}.svg",
                Self::sanitize_export_name_component(
                    alternate_reporter_seq_id,
                    "reporter_alternate"
                )
            ),
            report_md: "report.md".to_string(),
            result_json: "result.json".to_string(),
            commands_sh: "commands.sh".to_string(),
        }
    }

    fn variant_followup_context_viewport(
        &self,
        sequence_length_bp: usize,
    ) -> Option<(usize, usize)> {
        if sequence_length_bp == 0 {
            return None;
        }
        if let Ok(candidate) = self.variant_followup_recommended_candidate() {
            return Some((candidate.start_0based, candidate.length_bp.max(1)));
        }
        self.variant_followup_ui
            .cached_report
            .as_ref()
            .map(|report| {
                let flank = report.tfbs_focus_half_window_bp.max(200);
                let start = report.variant_start_0based.saturating_sub(flank);
                let end = report
                    .variant_end_0based_exclusive
                    .saturating_add(flank)
                    .min(sequence_length_bp);
                (
                    start.min(sequence_length_bp.saturating_sub(1)),
                    end.saturating_sub(start).max(1),
                )
            })
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

    fn start_variant_followup_score_tracks(&mut self) {
        let (seq_id, motifs, start_0based, end_0based_exclusive, score_kind, clip_negative) =
            match self.variant_followup_score_track_request() {
                Ok(value) => value,
                Err(err) => {
                    self.op_status = err;
                    return;
                }
            };
        self.start_tfbs_operation(
            Operation::SummarizeTfbsScoreTracks {
                target: SequenceScanTarget::SeqId {
                    seq_id,
                    span_start_0based: Some(start_0based),
                    span_end_0based_exclusive: Some(end_0based_exclusive),
                },
                motifs,
                score_kind,
                clip_negative,
                path: None,
            },
            TfbsTaskKind::VariantFollowupScoreTracks,
            "Promoter TF score tracks",
        );
    }

    fn summarize_variant_followup_score_tracks(&mut self) {
        self.start_variant_followup_score_tracks();
    }

    fn variant_followup_score_track_request(
        &self,
    ) -> Result<
        (
            String,
            Vec<String>,
            usize,
            usize,
            TfbsScoreTrackValueKind,
            bool,
        ),
        String,
    > {
        let seq_id = self.variant_followup_input_seq_id()?;
        let motifs =
            Self::promoter_design_parse_motif_tokens(&self.variant_followup_ui.score_track_motifs);
        if motifs.is_empty() {
            return Err(
                "Promoter design score tracks require at least one TF motif token".to_string(),
            );
        }
        let start_0based = Self::parse_optional_usize_text(
            &self.variant_followup_ui.score_track_start_0based,
            "score-track start bp",
        )?
        .unwrap_or(0);
        let end_0based_exclusive = Self::parse_optional_usize_text(
            &self.variant_followup_ui.score_track_end_0based_exclusive,
            "score-track end bp",
        )?
        .unwrap_or_else(|| {
            self.variant_followup_sequence_len(&seq_id)
                .unwrap_or(start_0based)
        });
        if start_0based >= end_0based_exclusive {
            return Err(format!(
                "Promoter design score tracks require start < end (got {}..{})",
                start_0based, end_0based_exclusive
            ));
        }
        Ok((
            seq_id,
            motifs,
            start_0based,
            end_0based_exclusive,
            self.variant_followup_ui.score_track_value_kind,
            self.variant_followup_ui.score_track_clip_negative,
        ))
    }

    fn export_variant_followup_score_tracks_svg(&mut self) {
        let (seq_id, motifs, start_0based, end_0based_exclusive, score_kind, clip_negative) =
            match self.variant_followup_score_track_request() {
                Ok(value) => value,
                Err(err) => {
                    self.op_status = err.clone();
                    self.op_error_popup = Some(err);
                    return;
                }
            };
        let default_name = format!(
            "{}_tfbs_score_tracks.svg",
            Self::sanitize_export_name_component(&seq_id, "promoter_design")
        );
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .save_file()
        else {
            self.op_status = "Promoter design TF score-track SVG export canceled".to_string();
            return;
        };
        let result =
            self.apply_operation_with_feedback_and_result(Operation::RenderTfbsScoreTracksSvg {
                target: SequenceScanTarget::SeqId {
                    seq_id,
                    span_start_0based: Some(start_0based),
                    span_end_0based_exclusive: Some(end_0based_exclusive),
                },
                motifs,
                score_kind,
                clip_negative,
                path: path.display().to_string(),
            });
        if let Some(report) = result.and_then(|row| row.tfbs_score_tracks) {
            self.variant_followup_ui.cached_score_tracks = Some(report);
        }
    }

    fn export_variant_followup_score_track_correlation_svg(&mut self) {
        let (seq_id, motifs, start_0based, end_0based_exclusive, score_kind, clip_negative) =
            match self.variant_followup_score_track_request() {
                Ok(value) => value,
                Err(err) => {
                    self.op_status = err.clone();
                    self.op_error_popup = Some(err);
                    return;
                }
            };
        let default_name = format!(
            "{}_tfbs_score_track_correlation.svg",
            Self::sanitize_export_name_component(&seq_id, "promoter_design")
        );
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .save_file()
        else {
            self.op_status = "Promoter design TF correlation SVG export canceled".to_string();
            return;
        };
        let result = self.apply_operation_with_feedback_and_result(
            Operation::RenderTfbsScoreTrackCorrelationSvg {
                seq_id,
                motifs,
                start_0based,
                end_0based_exclusive,
                score_kind,
                correlation_metric: self.variant_followup_ui.score_track_correlation_metric,
                correlation_signal_source: self
                    .variant_followup_ui
                    .score_track_correlation_signal_source,
                clip_negative,
                path: path.display().to_string(),
            },
        );
        if let Some(report) = result.and_then(|row| row.tfbs_score_tracks) {
            self.variant_followup_ui.cached_score_tracks = Some(report);
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

    fn run_variant_followup_promoter_context(
        &mut self,
        path: Option<String>,
    ) -> Option<VariantPromoterContextReport> {
        let input = match self.variant_followup_input_seq_id() {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return None;
            }
        };
        let promoter_upstream_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.promoter_upstream_bp,
            "promoter upstream bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return None;
            }
        };
        let promoter_downstream_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.promoter_downstream_bp,
            "promoter downstream bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return None;
            }
        };
        let tfbs_focus_half_window_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.tfbs_focus_half_window_bp,
            "TFBS focus half-window bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return None;
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
                path,
            },
        );
        if let Some(report) = result.and_then(|row| row.variant_promoter_context) {
            self.variant_followup_apply_report_defaults(&report);
            self.variant_followup_ui.cached_report = Some(report.clone());
            Some(report)
        } else {
            None
        }
    }

    fn summarize_variant_followup_promoter_context(&mut self) {
        let _ = self.run_variant_followup_promoter_context(None);
    }

    fn run_variant_followup_reporter_fragment_suggestion(
        &mut self,
        path: Option<String>,
    ) -> Option<PromoterReporterCandidateSet> {
        let input = match self.variant_followup_input_seq_id() {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return None;
            }
        };
        let retain_downstream_from_tss_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.retain_downstream_from_tss_bp,
            "retain downstream from TSS bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return None;
            }
        };
        let retain_upstream_beyond_variant_bp = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.retain_upstream_beyond_variant_bp,
            "retain upstream beyond variant bp",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return None;
            }
        };
        let max_candidates = match Self::parse_positive_usize_text(
            &self.variant_followup_ui.max_candidates,
            "max_candidates",
        ) {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return None;
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
                path,
            },
        );
        if let Some(candidates) = result.and_then(|row| row.promoter_reporter_candidates) {
            self.variant_followup_apply_candidate_defaults(&candidates);
            self.variant_followup_ui.cached_candidates = Some(candidates.clone());
            Some(candidates)
        } else {
            None
        }
    }

    fn suggest_variant_followup_reporter_fragments(&mut self) {
        let _ = self.run_variant_followup_reporter_fragment_suggestion(None);
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

    fn ensure_variant_followup_reporter_pair_preview_ids(
        &mut self,
    ) -> Result<(String, String), String> {
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
            return Err(
                "Matched promoter alleles are not available yet; extract and materialize them first"
                    .to_string(),
            );
        }
        if let Err(err) = self.ensure_variant_followup_reporter_backbone_loaded() {
            return Err(err);
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
        let reference_reporter_id = format!("{prefix}_reference");
        let alternate_reporter_id = format!("{prefix}_alternate");
        if self.variant_followup_sequence_exists(&reference_reporter_id)
            && self.variant_followup_sequence_exists(&alternate_reporter_id)
        {
            return Ok((reference_reporter_id, alternate_reporter_id));
        }
        let reference_assembly_prefix = format!("{prefix}_reference_assembly");
        let alternate_assembly_prefix = format!("{prefix}_alternate_assembly");

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
            return Err(
                "Reporter preview assembly did not yield a reference assembly id".to_string(),
            );
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
            return Err(
                "Reporter preview assembly did not yield a reference reporter id".to_string(),
            );
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
            return Err(
                "Reporter preview assembly did not yield an alternate assembly id".to_string(),
            );
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
            return Err(
                "Reporter preview assembly did not yield an alternate reporter id".to_string(),
            );
        };
        Ok((reference_preview_id, alternate_preview_id))
    }

    fn preview_variant_followup_reporter_pair(&mut self) {
        let (reference_preview_id, alternate_preview_id) =
            match self.ensure_variant_followup_reporter_pair_preview_ids() {
                Ok(value) => value,
                Err(err) => {
                    self.op_status = err.clone();
                    self.op_error_popup = Some(err);
                    return;
                }
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

    fn variant_followup_sequence_clone(&self, seq_id: &str) -> Result<DNAsequence, String> {
        if self.seq_id.as_deref() == Some(seq_id) {
            return Ok(self
                .dna
                .read()
                .map_err(|_| format!("Could not read active sequence '{seq_id}'"))?
                .clone());
        }
        let engine = self
            .engine
            .as_ref()
            .ok_or_else(|| "No engine attached".to_string())?;
        engine
            .read()
            .map_err(|_| "Engine lock poisoned".to_string())?
            .state()
            .sequences
            .get(seq_id)
            .cloned()
            .ok_or_else(|| format!("Sequence '{seq_id}' not found in engine state"))
    }

    fn variant_followup_export_display_settings(&self) -> DisplaySettings {
        let mut display = self
            .engine
            .as_ref()
            .and_then(|engine| {
                engine
                    .read()
                    .ok()
                    .map(|guard| guard.state().display.clone())
            })
            .unwrap_or_default();
        display.show_construct_reasoning_overlay = false;
        display.show_restriction_enzymes = false;
        display.show_gc_contents = false;
        display.show_open_reading_frames = false;
        display.show_methylation_sites = false;
        display
    }

    fn write_variant_followup_promoter_context_svg(
        &self,
        seq_id: &str,
        path: &Path,
    ) -> Result<(), String> {
        let dna = self.variant_followup_sequence_clone(seq_id)?;
        let mut display = self.variant_followup_export_display_settings();
        display.show_cds_features = false;
        display.show_gene_features = true;
        display.show_mrna_features = true;
        display.show_tfbs = true;
        if let Some((start_bp, span_bp)) = self.variant_followup_context_viewport(dna.len()) {
            display.linear_view_start_bp = start_bp;
            display.linear_view_span_bp = span_bp;
        }
        let svg = export_linear_svg(&dna, &display);
        fs::write(path, svg).map_err(|e| {
            format!(
                "Could not write promoter-context SVG '{}': {e}",
                path.display()
            )
        })
    }

    fn write_variant_followup_circular_svg(&self, seq_id: &str, path: &Path) -> Result<(), String> {
        let dna = self.variant_followup_sequence_clone(seq_id)?;
        let mut display = self.variant_followup_export_display_settings();
        display.show_tfbs = false;
        let svg = export_circular_svg(&dna, &display);
        fs::write(path, svg)
            .map_err(|e| format!("Could not write reporter SVG '{}': {e}", path.display()))
    }

    pub(super) fn build_variant_followup_handoff_result_json(
        &self,
        artifacts: &VariantFollowupBundleArtifacts,
        report: &VariantPromoterContextReport,
        candidates: &PromoterReporterCandidateSet,
        recommended_candidate: &crate::engine::PromoterReporterFragmentCandidate,
        reference_reporter_seq_id: &str,
        alternate_reporter_seq_id: &str,
    ) -> serde_json::Value {
        let genome_anchor = report.genome_anchor.as_ref();
        let variant_position_1based = genome_anchor.map(|anchor| {
            anchor
                .start_1based
                .saturating_add(report.variant_start_0based)
        });
        let variant_assembly = genome_anchor
            .map(|anchor| anchor.genome_id.clone())
            .unwrap_or_else(|| "unknown".to_string());
        serde_json::json!({
            "schema": "gentle.clawbio_handoff.v2",
            "id": artifacts.bundle_id.as_str(),
            "source_alert": {
                "drug": "warfarin",
                "gene": report.chosen_gene_label.clone().unwrap_or_else(|| "VKORC1".to_string()),
                "variant": {
                    "rs_id": report.variant_label.clone(),
                    "assembly": variant_assembly,
                    "chromosome": genome_anchor.map(|anchor| anchor.chromosome.clone()),
                    "position_1based": variant_position_1based
                }
            },
            "separation_of_concerns": {
                "clawbio": "Interpret the pharmacogenomic alert and motivate a regulatory follow-up around VKORC1/rs9923231.",
                "gentle": "Deterministically retrieve the locus, derive promoter context, suggest a reporter fragment, materialize matched alleles, load a mammalian reporter backbone, preview the constructs, and export reviewable artifacts."
            },
            "design": {
                "prepared_genome": genome_anchor.map(|anchor| anchor.genome_id.clone()),
                "context_sequence_id": report.seq_id.clone(),
                "promoter_window_defaults": {
                    "upstream_bp": report.promoter_upstream_bp,
                    "downstream_bp": report.promoter_downstream_bp,
                    "collapse_mode": "transcript"
                },
                "promoter_context_record": artifacts.variant_promoter_context_json.as_str(),
                "promoter_candidate_record": artifacts.promoter_reporter_candidates_json.as_str(),
                "fragment": {
                    "sequence_id": self.variant_followup_ui.fragment_output_id.trim(),
                    "engine_interval": {
                        "from": recommended_candidate.start_0based,
                        "to": recommended_candidate.end_0based_exclusive
                    },
                    "notes": [
                        format!(
                            "Selected from the recommended promoter-reporter candidate '{}'.",
                            recommended_candidate.candidate_id
                        ),
                        "Keeps the variant inside the insert rather than at the edge.",
                        format!(
                            "Retains {} bp downstream of the selected TSS and {} bp beyond the SNP on the biologically upstream side.",
                            recommended_candidate.retain_downstream_from_tss_bp,
                            recommended_candidate.retain_upstream_beyond_variant_bp
                        )
                    ]
                },
                "allele_specific_inserts": {
                    "reference": self.variant_followup_ui.reference_output_id.trim(),
                    "alternate": self.variant_followup_ui.alternate_output_id.trim()
                },
                "backbone": {
                    "sequence_id": self.variant_followup_ui.reporter_backbone_seq_id.trim(),
                    "path": self.variant_followup_ui.reporter_backbone_path.trim(),
                    "assay_role": "Pinned local promoterless mammalian luciferase reporter backbone for transient-transfection planning.",
                    "why_it_fits": [
                        "Keeps the assay focused on promoter output in human cells.",
                        "Supports matched upstream-of-luciferase promoter insert comparison.",
                        "Avoids making live GenBank retrieval part of the canonical tutorial path."
                    ]
                },
                "construct_previews": {
                    "status": "exported_from_variant_followup_expert",
                    "reference": {
                        "sequence_id": reference_reporter_seq_id,
                        "svg_path": artifacts.reference_reporter_svg.as_str()
                    },
                    "alternate": {
                        "sequence_id": alternate_reporter_seq_id,
                        "svg_path": artifacts.alternate_reporter_svg.as_str()
                    }
                },
                "promoter_context_svg": artifacts.promoter_context_svg.as_str(),
                "suggested_assay_family": if report.suggested_assay_ids.is_empty() {
                    candidates.suggested_assay_ids.first().cloned()
                } else {
                    report.suggested_assay_ids.first().cloned()
                }
            },
            "explicit_assumptions": [
                "This is a human-cell regulatory assay planning story, not a bacterial-expression story.",
                "The handoff ends at a reproducible construct-design point.",
                "Drug treatment is a later assay condition, not the first claim of this workflow.",
                "Adenoviral delivery is reserved as a later escalation only if transfection efficiency becomes the bottleneck."
            ],
            "artifacts": {
                "report_md": artifacts.report_md.as_str(),
                "result_json": artifacts.result_json.as_str(),
                "commands_sh": artifacts.commands_sh.as_str(),
                "promoter_context_svg": artifacts.promoter_context_svg.as_str(),
                "reference_reporter_svg": artifacts.reference_reporter_svg.as_str(),
                "alternate_reporter_svg": artifacts.alternate_reporter_svg.as_str()
            },
            "next_actions": [
                "Build the reference and alternate promoter fragments with identical boundaries.",
                "Keep the mammalian reporter backbone constant between the two constructs.",
                "Verify insert orientation and junction integrity before comparing reporter output.",
                "Choose one human cell model and normalization strategy for the later assay."
            ]
        })
    }

    fn build_variant_followup_handoff_report(
        &self,
        artifacts: &VariantFollowupBundleArtifacts,
        report: &VariantPromoterContextReport,
        candidates: &PromoterReporterCandidateSet,
        recommended_candidate: &crate::engine::PromoterReporterFragmentCandidate,
        reference_reporter_seq_id: &str,
        alternate_reporter_seq_id: &str,
    ) -> String {
        let chosen_gene = report
            .chosen_gene_label
            .clone()
            .or_else(|| candidates.chosen_gene_label.clone())
            .unwrap_or_else(|| "VKORC1".to_string());
        let chosen_transcript = report
            .chosen_transcript_label
            .clone()
            .or_else(|| report.chosen_transcript_id.clone())
            .or_else(|| candidates.chosen_transcript_label.clone())
            .or_else(|| candidates.chosen_transcript_id.clone())
            .unwrap_or_else(|| "<auto>".to_string());
        format!(
            "# Promoter design handoff\n\n## Interpretation coming from ClawBio\n\nClawBio interprets a pharmacogenomic alert around warfarin and `{variant}` and proposes a regulatory follow-up.\n\nThe narrow question handed to GENtle is:\n\n> Is `{variant}` promoter-proximal enough around `{gene}` that it is worth building matched promoter-reporter constructs for follow-up in human cells?\n\nClawBio is therefore responsible for:\n\n- the pharmacogenomic interpretation\n- the motivation for a regulatory assay\n- the choice to follow up in a human-cell luciferase reporter format\n\n## Deterministic build/render work done in GENtle\n\nGENtle is responsible for:\n\n- deriving transcript-TSS-centered promoter windows\n- summarizing promoter overlap and TSS distance\n- ranking transcript-aware reporter fragments\n- materializing matched reference and alternate inserts\n- loading a pinned mammalian promoterless luciferase backbone\n- previewing reference and alternate reporter constructs\n- exporting reviewable artifacts\n\n## Chosen baseline design\n\n- source context sequence id: `{source_seq}`\n- variant: `{variant}`\n- chosen gene: `{gene}`\n- chosen transcript: `{transcript}`\n- promoter overlap: `{promoter_overlap}`\n- signed TSS distance: `{signed_distance}`\n- promoter window defaults:\n  - upstream = `{promoter_upstream} bp`\n  - downstream = `{promoter_downstream} bp`\n- recommended promoter fragment interval:\n  - local `{fragment_start}..{fragment_end}`\n  - extracted id = `{fragment_id}`\n- matched inserts:\n  - `{reference_insert}`\n  - `{alternate_insert}`\n- pinned local backbone:\n  - `{backbone_id}`\n- construct previews:\n  - `{reference_reporter}`\n  - `{alternate_reporter}`\n\n## Why this backbone is appropriate for human-cell work\n\n- it is a promoterless mammalian luciferase reporter architecture\n- it keeps the readout focused on promoter output rather than bacterial protein-expression logic\n- it matches transient transfection planning in human cells\n- it gives the handoff one pinned local backbone so the canonical path does not depend on live GenBank retrieval\n\n## Important assumptions\n\n- this is a human-cell regulatory assay story\n- this bundle stops at a reproducible design/handoff point\n- it does **not** claim wet-lab validation\n- it does **not** claim that the construct alone proves warfarin response\n- adenoviral delivery is deferred as a later escalation only if transfection efficiency becomes the bottleneck\n\n## Artifacts\n\n- promoter-context JSON: `{promoter_context_json}`\n- promoter-candidate JSON: `{promoter_candidates_json}`\n- promoter-context SVG: `{promoter_context_svg}`\n- reference construct SVG: `{reference_svg}`\n- alternate construct SVG: `{alternate_svg}`\n- commands: `{commands_sh}`\n- structured summary: `{result_json}`\n\n## Bench-facing next actions\n\n1. Build the reference and alternate inserts with identical boundaries.\n2. Keep the mammalian reporter backbone constant between alleles.\n3. Verify insert orientation and junction integrity.\n4. Choose one human cell model and a normalization strategy for the later assay.\n5. Treat warfarin exposure as a later experimental condition layered on top of the finished promoter-reporter pair.\n",
            variant = report.variant_label.as_str(),
            gene = chosen_gene,
            source_seq = report.seq_id.as_str(),
            transcript = chosen_transcript,
            promoter_overlap = if report.promoter_overlap { "yes" } else { "no" },
            signed_distance = report
                .signed_tss_distance_bp
                .map(|value| format!("{value} bp"))
                .unwrap_or_else(|| "n/a".to_string()),
            promoter_upstream = report.promoter_upstream_bp,
            promoter_downstream = report.promoter_downstream_bp,
            fragment_start = recommended_candidate.start_0based,
            fragment_end = recommended_candidate.end_0based_exclusive,
            fragment_id = self.variant_followup_ui.fragment_output_id.trim(),
            reference_insert = self.variant_followup_ui.reference_output_id.trim(),
            alternate_insert = self.variant_followup_ui.alternate_output_id.trim(),
            backbone_id = self.variant_followup_ui.reporter_backbone_seq_id.trim(),
            reference_reporter = reference_reporter_seq_id,
            alternate_reporter = alternate_reporter_seq_id,
            promoter_context_json = artifacts.variant_promoter_context_json.as_str(),
            promoter_candidates_json = artifacts.promoter_reporter_candidates_json.as_str(),
            promoter_context_svg = artifacts.promoter_context_svg.as_str(),
            reference_svg = artifacts.reference_reporter_svg.as_str(),
            alternate_svg = artifacts.alternate_reporter_svg.as_str(),
            commands_sh = artifacts.commands_sh.as_str(),
            result_json = artifacts.result_json.as_str(),
        )
    }

    pub(super) fn build_variant_followup_handoff_commands(
        &self,
        artifacts: &VariantFollowupBundleArtifacts,
        report: &VariantPromoterContextReport,
        recommended_candidate: &crate::engine::PromoterReporterFragmentCandidate,
    ) -> String {
        let variant = report.variant_label.as_str();
        let source_seq_id = report.seq_id.as_str();
        let prepared_genome = report
            .genome_anchor
            .as_ref()
            .map(|anchor| anchor.genome_id.as_str())
            .unwrap_or("Human GRCh38 Ensembl 116");
        let flank_bp = report.variant_start_0based.max(
            report
                .sequence_length_bp
                .saturating_sub(report.variant_end_0based_exclusive),
        );
        let prefix = self.variant_followup_ui.reporter_output_prefix.trim();
        let prefix = if prefix.is_empty() {
            Self::variant_followup_suggested_token(variant)
        } else {
            Self::sanitize_export_name_component(prefix, "reporter_pair")
        };
        let reference_reporter_id = format!("{prefix}_reference");
        let alternate_reporter_id = format!("{prefix}_alternate");
        let reference_assembly_id = format!("{prefix}_reference_assembly_1");
        let alternate_assembly_id = format!("{prefix}_alternate_assembly_1");
        format!(
            "#!/usr/bin/env bash\nset -euo pipefail\n\nSTATE=\"${{STATE:-/tmp/{bundle_id}.state.json}}\"\nBUNDLE_DIR=\"${{BUNDLE_DIR:-$(pwd)}}\"\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  genomes status \"{prepared_genome}\" \\\n  --catalog assets/genomes.json \\\n  --cache-dir data/genomes\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  op '{{\"FetchDbSnpRegion\":{{\"rs_id\":\"{variant}\",\"genome_id\":\"{prepared_genome}\",\"flank_bp\":{flank_bp},\"output_id\":\"{source_seq_id}\",\"annotation_scope\":\"full\",\"catalog_path\":\"assets/genomes.json\",\"cache_dir\":\"data/genomes\"}}}}' \\\n  --confirm\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  variant annotate-promoters {source_seq_id} \\\n  --gene-label {gene_label} \\\n  --upstream-bp {promoter_upstream_bp} \\\n  --downstream-bp {promoter_downstream_bp}\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  variant promoter-context {source_seq_id} \\\n  --variant {variant} \\\n  --gene-label {gene_label} \\\n  --path \"$BUNDLE_DIR/{promoter_context_json}\"\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  variant reporter-fragments {source_seq_id} \\\n  --variant {variant} \\\n  --gene-label {gene_label} \\\n  --retain-downstream-from-tss-bp {retain_downstream_from_tss_bp} \\\n  --retain-upstream-beyond-variant-bp {retain_upstream_beyond_variant_bp} \\\n  --path \"$BUNDLE_DIR/{promoter_candidates_json}\"\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  op '{{\"ExtractRegion\":{{\"input\":\"{source_seq_id}\",\"from\":{fragment_start},\"to\":{fragment_end},\"output_id\":\"{fragment_id}\"}}}}' \\\n  --confirm\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  variant materialize-allele {fragment_id} \\\n  --variant {variant} \\\n  --allele reference \\\n  --output-id {reference_insert}\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  variant materialize-allele {fragment_id} \\\n  --variant {variant} \\\n  --allele alternate \\\n  --output-id {alternate_insert}\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  op '{{\"LoadFile\":{{\"path\":\"{backbone_path}\",\"as_id\":\"{backbone_id}\"}}}}' \\\n  --confirm\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  op '{{\"Ligation\":{{\"inputs\":[\"{reference_insert}\",\"{backbone_id}\"],\"circularize_if_possible\":false,\"protocol\":\"Blunt\",\"output_prefix\":\"{prefix}_reference_assembly\",\"unique\":false}}}}' \\\n  --confirm\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  op '{{\"Branch\":{{\"input\":\"{reference_assembly_id}\",\"output_id\":\"{reference_reporter_id}\"}}}}' \\\n  --confirm\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  op '{{\"Ligation\":{{\"inputs\":[\"{alternate_insert}\",\"{backbone_id}\"],\"circularize_if_possible\":false,\"protocol\":\"Blunt\",\"output_prefix\":\"{prefix}_alternate_assembly\",\"unique\":false}}}}' \\\n  --confirm\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  op '{{\"Branch\":{{\"input\":\"{alternate_assembly_id}\",\"output_id\":\"{alternate_reporter_id}\"}}}}' \\\n  --confirm\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  op '{{\"SetLinearViewport\":{{\"start_bp\":{viewport_start},\"span_bp\":{viewport_span}}}}}' \\\n  --confirm\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  render-svg {source_seq_id} linear \"$BUNDLE_DIR/{promoter_context_svg}\"\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  render-svg {reference_reporter_id} circular \"$BUNDLE_DIR/{reference_svg}\"\n\ncargo run --quiet --bin gentle_cli -- \\\n  --state \"$STATE\" \\\n  render-svg {alternate_reporter_id} circular \"$BUNDLE_DIR/{alternate_svg}\"\n",
            bundle_id = artifacts.bundle_id.as_str(),
            prepared_genome = prepared_genome,
            variant = variant,
            flank_bp = flank_bp,
            source_seq_id = source_seq_id,
            gene_label = Self::variant_followup_optional_text(&self.variant_followup_ui.gene_label)
                .unwrap_or_else(|| "VKORC1".to_string()),
            promoter_upstream_bp = report.promoter_upstream_bp,
            promoter_downstream_bp = report.promoter_downstream_bp,
            promoter_context_json = artifacts.variant_promoter_context_json.as_str(),
            retain_downstream_from_tss_bp = recommended_candidate.retain_downstream_from_tss_bp,
            retain_upstream_beyond_variant_bp =
                recommended_candidate.retain_upstream_beyond_variant_bp,
            promoter_candidates_json = artifacts.promoter_reporter_candidates_json.as_str(),
            fragment_start = recommended_candidate.start_0based,
            fragment_end = recommended_candidate.end_0based_exclusive,
            fragment_id = self.variant_followup_ui.fragment_output_id.trim(),
            reference_insert = self.variant_followup_ui.reference_output_id.trim(),
            alternate_insert = self.variant_followup_ui.alternate_output_id.trim(),
            backbone_path = self.variant_followup_ui.reporter_backbone_path.trim(),
            backbone_id = self.variant_followup_ui.reporter_backbone_seq_id.trim(),
            prefix = prefix,
            reference_assembly_id = reference_assembly_id,
            reference_reporter_id = reference_reporter_id,
            alternate_assembly_id = alternate_assembly_id,
            alternate_reporter_id = alternate_reporter_id,
            viewport_start = self
                .variant_followup_context_viewport(report.sequence_length_bp)
                .map(|(start, _)| start)
                .unwrap_or(recommended_candidate.start_0based),
            viewport_span = self
                .variant_followup_context_viewport(report.sequence_length_bp)
                .map(|(_, span)| span)
                .unwrap_or(recommended_candidate.length_bp.max(1)),
            promoter_context_svg = artifacts.promoter_context_svg.as_str(),
            reference_svg = artifacts.reference_reporter_svg.as_str(),
            alternate_svg = artifacts.alternate_reporter_svg.as_str(),
        )
    }

    fn export_variant_followup_handoff_bundle(&mut self) {
        let Some(parent_dir) = rfd::FileDialog::new().pick_folder() else {
            self.op_status = "Promoter design bundle export canceled".to_string();
            return;
        };
        let report = match self.run_variant_followup_promoter_context(None) {
            Some(value) => value,
            None => {
                if self.op_status.is_empty() {
                    self.op_status =
                        "Could not export handoff bundle because promoter context failed"
                            .to_string();
                }
                return;
            }
        };
        let candidates = match self.run_variant_followup_reporter_fragment_suggestion(None) {
            Some(value) => value,
            None => {
                if self.op_status.is_empty() {
                    self.op_status =
                        "Could not export handoff bundle because reporter-fragment suggestion failed"
                            .to_string();
                }
                return;
            }
        };
        let recommended_candidate = match self.variant_followup_recommended_candidate() {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err.clone();
                self.op_error_popup = Some(err);
                return;
            }
        };
        self.extract_variant_followup_recommended_fragment();
        if self
            .variant_followup_ui
            .fragment_output_id
            .trim()
            .is_empty()
            || !self.variant_followup_sequence_exists(
                self.variant_followup_ui.fragment_output_id.trim(),
            )
        {
            return;
        }
        self.materialize_variant_followup_alleles();
        let (reference_reporter_seq_id, alternate_reporter_seq_id) =
            match self.ensure_variant_followup_reporter_pair_preview_ids() {
                Ok(value) => value,
                Err(err) => {
                    self.op_status = err.clone();
                    self.op_error_popup = Some(err);
                    return;
                }
            };
        let artifacts = self.variant_followup_bundle_artifacts(
            &reference_reporter_seq_id,
            &alternate_reporter_seq_id,
        );
        let bundle_dir = parent_dir.join(&artifacts.bundle_id);
        if let Err(err) = fs::create_dir_all(&bundle_dir) {
            self.op_status = format!(
                "Could not create Promoter design bundle directory '{}': {err}",
                bundle_dir.display()
            );
            return;
        }

        let promoter_context_json_path = bundle_dir.join(&artifacts.variant_promoter_context_json);
        if self
            .run_variant_followup_promoter_context(Some(
                promoter_context_json_path.display().to_string(),
            ))
            .is_none()
        {
            return;
        }
        let promoter_candidates_json_path =
            bundle_dir.join(&artifacts.promoter_reporter_candidates_json);
        if self
            .run_variant_followup_reporter_fragment_suggestion(Some(
                promoter_candidates_json_path.display().to_string(),
            ))
            .is_none()
        {
            return;
        }

        let promoter_context_svg_path = bundle_dir.join(&artifacts.promoter_context_svg);
        if let Err(err) = self
            .write_variant_followup_promoter_context_svg(&report.seq_id, &promoter_context_svg_path)
        {
            self.op_status = err.clone();
            self.op_error_popup = Some(err);
            return;
        }
        let reference_svg_path = bundle_dir.join(&artifacts.reference_reporter_svg);
        if let Err(err) = self
            .write_variant_followup_circular_svg(&reference_reporter_seq_id, &reference_svg_path)
        {
            self.op_status = err.clone();
            self.op_error_popup = Some(err);
            return;
        }
        let alternate_svg_path = bundle_dir.join(&artifacts.alternate_reporter_svg);
        if let Err(err) = self
            .write_variant_followup_circular_svg(&alternate_reporter_seq_id, &alternate_svg_path)
        {
            self.op_status = err.clone();
            self.op_error_popup = Some(err);
            return;
        }

        let result_json = self.build_variant_followup_handoff_result_json(
            &artifacts,
            &report,
            &candidates,
            &recommended_candidate,
            &reference_reporter_seq_id,
            &alternate_reporter_seq_id,
        );
        let report_md = self.build_variant_followup_handoff_report(
            &artifacts,
            &report,
            &candidates,
            &recommended_candidate,
            &reference_reporter_seq_id,
            &alternate_reporter_seq_id,
        );
        let commands_sh = self.build_variant_followup_handoff_commands(
            &artifacts,
            &report,
            &recommended_candidate,
        );
        let result_json_path = bundle_dir.join(&artifacts.result_json);
        let report_md_path = bundle_dir.join(&artifacts.report_md);
        let commands_sh_path = bundle_dir.join(&artifacts.commands_sh);
        let result_json_text = match serde_json::to_string_pretty(&result_json) {
            Ok(text) => text,
            Err(err) => {
                self.op_status = format!("Could not serialize Promoter design result.json: {err}");
                return;
            }
        };
        for (path, text, label) in [
            (&result_json_path, result_json_text.as_str(), "result.json"),
            (&report_md_path, report_md.as_str(), "report.md"),
            (&commands_sh_path, commands_sh.as_str(), "commands.sh"),
        ] {
            if let Err(err) = fs::write(path, text) {
                self.op_status = format!(
                    "Could not write Promoter design {label} '{}': {err}",
                    path.display()
                );
                return;
            }
        }

        self.last_created_seq_ids = vec![
            reference_reporter_seq_id.clone(),
            alternate_reporter_seq_id.clone(),
        ];
        self.export_pool_inputs_text = self.last_created_seq_ids.join(", ");
        self.op_status = format!(
            "Exported Promoter design bundle to '{}'",
            bundle_dir.display()
        );
        self.op_error_popup = None;
    }

    fn render_variant_followup_report_summary(&mut self, ui: &mut egui::Ui) {
        if !self.variant_followup_has_variant_seed() {
            ui.small(
                egui::RichText::new(
                    "No variant is selected for this Promoter design session yet. Open the window from a variation feature to capture promoter-context reasoning and allele-specific reporter planning.",
                )
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
            return;
        }
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
        if !self.variant_followup_has_variant_seed() {
            ui.small(
                egui::RichText::new(
                    "Reporter-fragment ranking is variant-specific. Open Promoter design from a variation feature when you want GENtle to propose allele-paired reporter inserts.",
                )
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
            return;
        }
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

    fn promoter_design_plot_markers(&self) -> Vec<(usize, &'static str, egui::Color32)> {
        let mut markers = vec![];
        if let Some(report) = self.variant_followup_ui.cached_report.as_ref() {
            markers.push((
                report.variant_start_0based,
                "variant",
                egui::Color32::from_rgb(245, 158, 11),
            ));
        } else if self.variant_followup_has_variant_seed()
            && let (Some(source_seq_id), Some(feature_id)) = (
                Some(self.variant_followup_ui.source_seq_id.as_str()),
                self.variant_followup_ui.source_feature_id,
            )
            && let Some(feature) = self.variant_followup_feature_clone(source_seq_id, feature_id)
        {
            let mut ranges = vec![];
            collect_location_ranges_usize(&feature.location, &mut ranges);
            if let Some((start, _)) = ranges.into_iter().min_by_key(|(start, _)| *start) {
                markers.push((start, "variant", egui::Color32::from_rgb(245, 158, 11)));
            }
        }
        markers
    }

    fn promoter_design_track_normalization_summary(
        score_kind: TfbsScoreTrackValueKind,
        normalization: &crate::engine::TfbsScoreTrackNormalizationReference,
    ) -> String {
        match score_kind {
            TfbsScoreTrackValueKind::LlrBackgroundQuantile
            | TfbsScoreTrackValueKind::LlrBackgroundTailLog10
            | TfbsScoreTrackValueKind::TrueLogOddsBackgroundQuantile
            | TfbsScoreTrackValueKind::TrueLogOddsBackgroundTailLog10 => format!(
                "theory max {:.2} | peak q {:.6} | -log10 tail {:.2}",
                normalization.theoretical_max_score,
                normalization.observed_peak_modeled_quantile,
                normalization.observed_peak_modeled_tail_log10,
            ),
            _ => format!(
                "p99 {:.2} | Δp99 {:+.2} | bg+ {:.1}%",
                normalization.p99_score,
                normalization.observed_peak_delta_from_p99,
                normalization.positive_fraction.max(0.0) * 100.0
            ),
        }
    }

    pub(super) fn promoter_design_track_label(tf_id: &str, tf_name: Option<&str>) -> String {
        let trimmed_name = tf_name.map(str::trim).unwrap_or_default();
        if trimmed_name.is_empty() {
            return tf_id.to_string();
        }
        if trimmed_name.eq_ignore_ascii_case(tf_id) {
            trimmed_name.to_string()
        } else {
            format!("{trimmed_name} ({tf_id})")
        }
    }

    pub(super) fn render_tfbs_task_progress_panel(
        ui: &mut egui::Ui,
        task: &TfbsTask,
        progress: Option<&TfbsProgress>,
        cancelable: bool,
    ) -> bool {
        let mut cancel_clicked = false;
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.add(egui::Spinner::new());
                ui.label(MainAreaDna::summarize_tfbs_task_status(task, progress));
            });
            if cancelable {
                ui.horizontal(|ui| {
                    let cancel_requested = task.cancel_requested.load(AtomicOrdering::Relaxed);
                    if cancel_requested {
                        ui.small(
                            egui::RichText::new(
                                "Cancel requested... waiting for TFBS worker to stop.",
                            )
                            .color(egui::Color32::from_rgb(220, 38, 38)),
                        );
                    } else if ui
                        .button(format!("Cancel {}", task.operation_label))
                        .on_hover_text("Request cancellation of the running TFBS computation.")
                        .clicked()
                    {
                        cancel_clicked = true;
                    }
                });
            }
            if let Some(progress) = progress {
                ui.add(
                    egui::ProgressBar::new((progress.total_percent / 100.0) as f32)
                        .show_percentage()
                        .text(format!(
                            "Overall: motif {}/{}",
                            progress.motif_index, progress.motif_count
                        )),
                );
                ui.add(
                    egui::ProgressBar::new(
                        (progress.stage_percent.unwrap_or(progress.motif_percent) / 100.0) as f32,
                    )
                    .show_percentage()
                    .text(MainAreaDna::tfbs_progress_detail_label(progress)),
                );
            }
        });
        cancel_clicked
    }

    fn promoter_design_track_peak_summary(
        track: &crate::engine::TfbsScoreTrackRow,
    ) -> Option<String> {
        if track.top_peaks.is_empty() {
            return None;
        }
        Some(format!(
            "top peaks: {}",
            track
                .top_peaks
                .iter()
                .map(|peak| format!(
                    "{}{} {:.2} (Δp99 {:+.2})",
                    peak.start_0based,
                    if peak.is_reverse { "r" } else { "f" },
                    peak.score,
                    peak.delta_from_p99
                ))
                .collect::<Vec<_>>()
                .join(", ")
        ))
    }

    fn render_variant_followup_score_track_correlation_summary(
        ui: &mut egui::Ui,
        report: &TfbsScoreTrackReport,
        correlation_metric: Option<&mut TfbsScoreTrackCorrelationMetric>,
        correlation_signal_source: Option<&mut TfbsScoreTrackCorrelationSignalSource>,
    ) {
        ui.add_space(8.0);
        ui.group(|ui| {
            ui.label(egui::RichText::new("TFBS track correlation").strong());
            let metric = match correlation_metric {
                Some(metric) => {
                    ui.horizontal(|ui| {
                        ui.small("metric");
                        egui::ComboBox::from_id_salt(
                            "promoter_design_score_track_correlation_metric",
                        )
                        .selected_text(metric.as_str())
                        .show_ui(ui, |ui| {
                            for choice in [
                                TfbsScoreTrackCorrelationMetric::Pearson,
                                TfbsScoreTrackCorrelationMetric::Spearman,
                            ] {
                                ui.selectable_value(metric, choice, choice.as_str());
                            }
                        });
                    });
                    *metric
                }
                None => TfbsScoreTrackCorrelationMetric::Pearson,
            };
            let signal_source = match correlation_signal_source {
                Some(signal_source) => {
                    ui.horizontal(|ui| {
                        ui.small("signal");
                        egui::ComboBox::from_id_salt(
                            "promoter_design_score_track_correlation_signal_source",
                        )
                        .selected_text(signal_source.as_str())
                        .show_ui(ui, |ui| {
                            for choice in [
                                TfbsScoreTrackCorrelationSignalSource::MaxStrands,
                                TfbsScoreTrackCorrelationSignalSource::ForwardOnly,
                                TfbsScoreTrackCorrelationSignalSource::ReverseOnly,
                            ] {
                                ui.selectable_value(
                                    signal_source,
                                    choice,
                                    choice.as_str(),
                                );
                            }
                        });
                    });
                    *signal_source
                }
                None => TfbsScoreTrackCorrelationSignalSource::MaxStrands,
            };
            let correlation = report
                .correlation_summaries
                .iter()
                .find(|summary| summary.signal_source == signal_source)
                .or_else(|| {
                    report.correlation_summary.as_ref().filter(|summary| {
                        summary.signal_source == signal_source
                            || signal_source == TfbsScoreTrackCorrelationSignalSource::MaxStrands
                    })
                });
            let Some(correlation) = correlation else {
                ui.small(
                    egui::RichText::new(
                        "No pairwise correlation summary is available for this signal source yet.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
                return;
            };
            if correlation.rows.is_empty() {
                ui.small(
                    egui::RichText::new(
                        "No pairwise correlation rows are available for the selected signal source.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
                return;
            }
            ui.small(
                egui::RichText::new(format!(
                    "Signal source: {} | smoothed view: {} {} bp | raw and smoothed {} are shown.",
                    correlation.signal_source.summary_label(),
                    correlation.smoothing_method,
                    correlation.smoothing_window_bp,
                    metric.as_str()
                ))
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
            let directional_rows = report
                .tracks
                .iter()
                .filter_map(|track| {
                    track.directional_summary.as_ref().map(|summary| (track, summary))
                })
                .collect::<Vec<_>>();
            if !directional_rows.is_empty() {
                ui.add_space(4.0);
                ui.small(
                    egui::RichText::new(
                        "Per-motif strand pairing shows how forward and reverse score tracks travel together inside the same motif row.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
                egui::Grid::new((
                    "variant_followup_tfbs_track_directional_grid",
                    report.seq_id.as_str(),
                    report.score_kind.as_str(),
                ))
                .num_columns(4)
                .spacing([12.0, 4.0])
                .show(ui, |ui| {
                    ui.small("Motif");
                    ui.small(metric.display_label(true));
                    ui.small(metric.display_label(false));
                    ui.small("Rev-Fwd peak offset");
                    ui.end_row();
                    for (track, summary) in directional_rows {
                        ui.small(Self::promoter_design_track_label(
                            &track.tf_id,
                            track.tf_name.as_deref(),
                        ));
                        ui.monospace(format!(
                            "{:+.3}",
                            match metric {
                                TfbsScoreTrackCorrelationMetric::Pearson => {
                                    summary.smoothed_pearson
                                }
                                TfbsScoreTrackCorrelationMetric::Spearman => {
                                    summary.smoothed_spearman
                                }
                            }
                        ));
                        ui.monospace(format!(
                            "{:+.3}",
                            match metric {
                                TfbsScoreTrackCorrelationMetric::Pearson => summary.raw_pearson,
                                TfbsScoreTrackCorrelationMetric::Spearman => summary.raw_spearman,
                            }
                        ));
                        ui.monospace(
                            summary
                                .signed_primary_peak_offset_bp
                                .map(|value| format!("{value:+} bp"))
                                .unwrap_or_else(|| "n/a".to_string()),
                        );
                        ui.end_row();
                    }
                });
            }
            ui.add_space(4.0);
            egui::Grid::new((
                "variant_followup_tfbs_track_correlation_grid",
                report.seq_id.as_str(),
                report.score_kind.as_str(),
            ))
            .num_columns(4)
            .spacing([12.0, 4.0])
            .show(ui, |ui| {
                let mut rows = correlation.rows.iter().collect::<Vec<_>>();
                rows.sort_by(|left, right| {
                    let (right_smoothed, right_raw) = match metric {
                        TfbsScoreTrackCorrelationMetric::Pearson => {
                            (right.smoothed_pearson.abs(), right.raw_pearson.abs())
                        }
                        TfbsScoreTrackCorrelationMetric::Spearman => {
                            (right.smoothed_spearman.abs(), right.raw_spearman.abs())
                        }
                    };
                    let (left_smoothed, left_raw) = match metric {
                        TfbsScoreTrackCorrelationMetric::Pearson => {
                            (left.smoothed_pearson.abs(), left.raw_pearson.abs())
                        }
                        TfbsScoreTrackCorrelationMetric::Spearman => {
                            (left.smoothed_spearman.abs(), left.raw_spearman.abs())
                        }
                    };
                    right_smoothed
                        .total_cmp(&left_smoothed)
                        .then(right_raw.total_cmp(&left_raw))
                        .then(left.left_tf_id.cmp(&right.left_tf_id))
                        .then(left.right_tf_id.cmp(&right.right_tf_id))
                });
                ui.small("Pair");
                ui.small(metric.display_label(true));
                ui.small(metric.display_label(false));
                ui.small("Peak offset");
                ui.end_row();
                for row in rows {
                    let left_label = Self::promoter_design_track_label(
                        &row.left_tf_id,
                        row.left_tf_name.as_deref(),
                    );
                    let right_label = Self::promoter_design_track_label(
                        &row.right_tf_id,
                        row.right_tf_name.as_deref(),
                    );
                    ui.small(format!("{left_label} <> {right_label}"));
                    ui.monospace(format!(
                        "{:+.3}",
                        match metric {
                            TfbsScoreTrackCorrelationMetric::Pearson => row.smoothed_pearson,
                            TfbsScoreTrackCorrelationMetric::Spearman => row.smoothed_spearman,
                        }
                    ));
                    ui.monospace(format!(
                        "{:+.3}",
                        match metric {
                            TfbsScoreTrackCorrelationMetric::Pearson => row.raw_pearson,
                            TfbsScoreTrackCorrelationMetric::Spearman => row.raw_spearman,
                        }
                    ));
                    ui.monospace(
                        row.signed_primary_peak_offset_bp
                            .map(|value| format!("{value:+} bp"))
                            .unwrap_or_else(|| "n/a".to_string()),
                    );
                    ui.end_row();
                }
            });
        });
    }

    fn paint_promoter_design_track_plot(
        ui: &mut egui::Ui,
        report: &TfbsScoreTrackReport,
        track: &crate::engine::TfbsScoreTrackRow,
        markers: &[(usize, &'static str, egui::Color32)],
    ) {
        let desired = egui::vec2(ui.available_width().max(420.0), 76.0);
        let (rect, _) = ui.allocate_exact_size(desired, egui::Sense::hover());
        let painter = ui.painter_at(rect);
        let plot_rect = rect.shrink2(egui::vec2(8.0, 8.0));
        painter.rect_stroke(
            plot_rect,
            4.0,
            egui::Stroke::new(1.0, egui::Color32::from_rgb(203, 213, 225)),
            egui::StrokeKind::Inside,
        );
        let baseline_y = plot_rect.bottom() - 10.0;
        painter.line_segment(
            [
                egui::pos2(plot_rect.left(), baseline_y),
                egui::pos2(plot_rect.right(), baseline_y),
            ],
            egui::Stroke::new(1.0, egui::Color32::from_rgb(148, 163, 184)),
        );
        if track.scored_window_count == 0 {
            painter.text(
                plot_rect.center(),
                egui::Align2::CENTER_CENTER,
                "motif longer than selected range",
                egui::FontId::proportional(12.0),
                egui::Color32::from_rgb(100, 116, 139),
            );
            return;
        }
        let score_max = report.global_max_score.max(1.0) as f32;
        let x_denom = track.scored_window_count.saturating_sub(1).max(1) as f32;
        let y_height = (baseline_y - (plot_rect.top() + 8.0)).max(8.0);
        let to_x = |idx: usize| plot_rect.left() + (idx as f32 / x_denom) * plot_rect.width();
        let to_y =
            |score: f64| baseline_y - ((score as f32 / score_max).clamp(0.0, 1.0) * y_height);

        for (marker_pos, label, color) in markers {
            if *marker_pos < report.view_start_0based
                || *marker_pos >= report.view_end_0based_exclusive
            {
                continue;
            }
            let span = report
                .view_end_0based_exclusive
                .saturating_sub(report.view_start_0based)
                .saturating_sub(1)
                .max(1) as f32;
            let x = plot_rect.left()
                + ((*marker_pos - report.view_start_0based) as f32 / span) * plot_rect.width();
            painter.line_segment(
                [
                    egui::pos2(x, plot_rect.top() + 4.0),
                    egui::pos2(x, baseline_y),
                ],
                egui::Stroke::new(1.0, *color),
            );
            painter.text(
                egui::pos2(x + 3.0, plot_rect.top() + 2.0),
                egui::Align2::LEFT_TOP,
                *label,
                egui::FontId::proportional(10.0),
                *color,
            );
        }

        let forward_points = track
            .forward_scores
            .iter()
            .enumerate()
            .map(|(idx, score)| egui::pos2(to_x(idx), to_y(*score)))
            .collect::<Vec<_>>();
        if forward_points.len() > 1 {
            painter.add(egui::Shape::line(
                forward_points,
                egui::Stroke::new(1.8, egui::Color32::from_rgb(14, 116, 144)),
            ));
        }
        let reverse_points = track
            .reverse_scores
            .iter()
            .enumerate()
            .map(|(idx, score)| egui::pos2(to_x(idx), to_y(*score)))
            .collect::<Vec<_>>();
        if reverse_points.len() > 1 {
            painter.add(egui::Shape::line(
                reverse_points,
                egui::Stroke::new(1.6, egui::Color32::from_rgb(180, 83, 9)),
            ));
        }
        painter.text(
            egui::pos2(plot_rect.right() - 4.0, plot_rect.top() + 2.0),
            egui::Align2::RIGHT_TOP,
            format!("0 .. {:.2}", report.global_max_score),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(71, 85, 105),
        );
    }

    fn promoter_design_overlay_track_meta(
        track: &crate::engine::TfbsScoreTrackOverlayTrack,
    ) -> String {
        let mut parts = vec![format!("{} interval(s)", track.interval_count)];
        if let Some(max_score) = track.max_score {
            parts.push(format!("max score {:.2}", max_score));
        }
        format!("{} | {}", track.source_kind, parts.join(" | "))
    }

    fn promoter_design_overlay_interval_alpha(
        interval: &crate::engine::TfbsScoreTrackOverlayInterval,
        track: &crate::engine::TfbsScoreTrackOverlayTrack,
    ) -> f32 {
        match (interval.score, track.max_score) {
            (Some(score), Some(max_score)) if max_score.is_finite() && max_score > 0.0 => {
                ((score / max_score).clamp(0.0, 1.0) as f32 * 0.5 + 0.25).clamp(0.25, 0.85)
            }
            _ => 0.6,
        }
    }

    fn paint_promoter_design_overlay_track_plot(
        ui: &mut egui::Ui,
        report: &TfbsScoreTrackReport,
        track: &crate::engine::TfbsScoreTrackOverlayTrack,
    ) {
        let desired = egui::vec2(ui.available_width().max(420.0), 34.0);
        let (rect, _) = ui.allocate_exact_size(desired, egui::Sense::hover());
        let painter = ui.painter_at(rect);
        let plot_rect = rect.shrink2(egui::vec2(8.0, 6.0));
        painter.rect_filled(plot_rect, 4.0, egui::Color32::WHITE);
        painter.rect_stroke(
            plot_rect,
            4.0,
            egui::Stroke::new(1.0, egui::Color32::from_rgb(221, 214, 254)),
            egui::StrokeKind::Inside,
        );
        let denom = report
            .view_end_0based_exclusive
            .saturating_sub(report.view_start_0based)
            .max(1) as f32;
        let lane_rect = plot_rect.shrink2(egui::vec2(6.0, 5.0));
        for interval in &track.intervals {
            let start_fraction = (interval
                .start_0based
                .saturating_sub(report.view_start_0based) as f32
                / denom)
                .clamp(0.0, 1.0);
            let end_fraction = (interval
                .end_0based_exclusive
                .saturating_sub(report.view_start_0based) as f32
                / denom)
                .clamp(0.0, 1.0);
            let x1 = lane_rect.left() + lane_rect.width() * start_fraction;
            let x2 = lane_rect.left() + lane_rect.width() * end_fraction;
            let width = (x2 - x1).max(2.0);
            let fill = egui::Color32::from_rgba_premultiplied(
                124,
                58,
                237,
                (255.0 * Self::promoter_design_overlay_interval_alpha(interval, track)) as u8,
            );
            let interval_rect = egui::Rect::from_min_size(
                egui::pos2(x1, lane_rect.top()),
                egui::vec2(width, lane_rect.height()),
            );
            painter.rect_filled(interval_rect, 2.0, fill);
            painter.rect_stroke(
                interval_rect,
                2.0,
                egui::Stroke::new(0.8, egui::Color32::from_rgb(91, 33, 182)),
                egui::StrokeKind::Inside,
            );
        }
    }

    pub(super) fn render_tfbs_score_track_summary_panel(
        ui: &mut egui::Ui,
        title: &str,
        empty_message: &str,
        report: Option<&TfbsScoreTrackReport>,
        markers: &[(usize, &'static str, egui::Color32)],
        correlation_metric: Option<&mut TfbsScoreTrackCorrelationMetric>,
        correlation_signal_source: Option<&mut TfbsScoreTrackCorrelationSignalSource>,
    ) {
        let Some(report) = report else {
            ui.small(
                egui::RichText::new(empty_message).color(egui::Color32::from_rgb(100, 116, 139)),
            );
            return;
        };
        let target_label = if report.target_label.trim().is_empty() {
            report.seq_id.as_str()
        } else {
            report.target_label.trim()
        };
        ui.group(|ui| {
            ui.label(egui::RichText::new(title).strong());
            ui.small(
                egui::RichText::new(format!(
                    "{} motif(s) across {}..{} | score={}{}",
                    report.tracks.len(),
                    report.view_start_0based,
                    report.view_end_0based_exclusive,
                    report.score_kind.as_str(),
                    if report.clip_negative && report.score_kind.supports_negative_values() {
                        " | negative values clipped to 0"
                    } else {
                        ""
                    }
                ))
                .color(egui::Color32::from_rgb(71, 85, 105)),
            );
            ui.small(
                egui::RichText::new(format!(
                    "target: {} '{}' | source length {} bp | topology {}",
                    report.target_kind,
                    target_label,
                    report.source_sequence_length_bp,
                    report.scan_topology.as_str()
                ))
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
            if let Some(normalization) = report
                .tracks
                .iter()
                .find_map(|track| track.normalization_reference.as_ref())
            {
                ui.small(
                    egui::RichText::new(format!(
                        "background normalization: {} {} bp deterministic random DNA",
                        normalization.background_model, normalization.random_sequence_length_bp
                    ))
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
            }
            for track in &report.tracks {
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.set_min_width(238.0);
                    ui.horizontal(|ui| {
                        ui.vertical(|ui| {
                            let label = Self::promoter_design_track_label(
                                &track.tf_id,
                                track.tf_name.as_deref(),
                            );
                            ui.horizontal_wrapped(|ui| {
                                ui.label(egui::RichText::new(label).strong());
                                if ui
                                    .small_button("PSSM…")
                                    .on_hover_text(
                                        "Open the JASPAR Expert focused on this motif/PSSM.",
                                    )
                                    .clicked()
                                {
                                    crate::app::request_open_jaspar_expert_for_motif_from_native_menu(
                                        &track.tf_id,
                                    );
                                }
                            });
                            ui.small(format!(
                                "{} bp motif | {} windows | max {:.2}{}",
                                track.motif_length_bp,
                                track.scored_window_count,
                                track.max_score,
                                track
                                    .max_position_0based
                                    .map(|pos| format!(" @ {}", pos))
                                    .unwrap_or_default()
                            ));
                            if let Some(normalization) = track.normalization_reference.as_ref() {
                                ui.small(
                                    egui::RichText::new(
                                        Self::promoter_design_track_normalization_summary(
                                            report.score_kind,
                                            normalization,
                                        ),
                                    )
                                    .color(egui::Color32::from_rgb(100, 116, 139)),
                                );
                            }
                            if let Some(peaks_text) =
                                Self::promoter_design_track_peak_summary(track)
                            {
                                ui.small(
                                    egui::RichText::new(peaks_text)
                                        .color(egui::Color32::from_rgb(71, 85, 105)),
                                );
                            }
                            ui.small(
                                egui::RichText::new("forward / reverse")
                                    .color(egui::Color32::from_rgb(71, 85, 105)),
                            );
                        });
                    });
                    ui.vertical(|ui| {
                        Self::paint_promoter_design_track_plot(ui, report, track, &markers);
                    });
                });
            }
            if !report.overlay_tracks.is_empty() {
                ui.add_space(6.0);
                ui.label(
                    egui::RichText::new("Imported BED interval tracks")
                        .color(egui::Color32::from_rgb(88, 28, 135))
                        .strong(),
                );
                ui.small(
                    egui::RichText::new(
                        "Imported BED/CUT&RUN peaks on this sequence are shown as extra interval lanes under the motif score tracks.",
                    )
                    .color(egui::Color32::from_rgb(107, 114, 128)),
                );
                for overlay_track in &report.overlay_tracks {
                    ui.add_space(4.0);
                    ui.horizontal(|ui| {
                        ui.set_min_width(238.0);
                        ui.vertical(|ui| {
                            ui.label(
                                egui::RichText::new(&overlay_track.display_label)
                                    .color(egui::Color32::from_rgb(88, 28, 135))
                                    .strong(),
                            );
                            ui.small(
                                egui::RichText::new(
                                    Self::promoter_design_overlay_track_meta(overlay_track),
                                )
                                .color(egui::Color32::from_rgb(107, 114, 128)),
                            );
                        });
                        ui.vertical(|ui| {
                            Self::paint_promoter_design_overlay_track_plot(
                                ui,
                                report,
                                overlay_track,
                            );
                        });
                    });
                }
            }
            Self::render_variant_followup_score_track_correlation_summary(
                ui,
                report,
                correlation_metric,
                correlation_signal_source,
            );
        });
    }

    fn render_variant_followup_score_track_summary(&mut self, ui: &mut egui::Ui) {
        let markers = self.promoter_design_plot_markers();
        let mut correlation_metric = self.variant_followup_ui.score_track_correlation_metric;
        let mut correlation_signal_source = self
            .variant_followup_ui
            .score_track_correlation_signal_source;
        if let Some(task) = self
            .tfbs_task
            .as_ref()
            .filter(|task| matches!(task.task_kind, TfbsTaskKind::VariantFollowupScoreTracks))
        {
            if Self::render_tfbs_task_progress_panel(ui, task, self.tfbs_progress.as_ref(), true) {
                task.cancel_requested.store(true, AtomicOrdering::Relaxed);
                self.op_status = format!("Cancel requested for {}", task.operation_label);
            }
            ui.add_space(6.0);
        }
        Self::render_tfbs_score_track_summary_panel(
            ui,
            "Promoter TF score tracks",
            "No TF score tracks cached yet. Run 'Show TF score tracks' to inspect positive-only motif scoring across the selected promoter span.",
            self.variant_followup_ui.cached_score_tracks.as_ref(),
            &markers,
            Some(&mut correlation_metric),
            Some(&mut correlation_signal_source),
        );
        self.variant_followup_ui.score_track_correlation_metric = correlation_metric;
        self.variant_followup_ui
            .score_track_correlation_signal_source = correlation_signal_source;
    }

    fn render_variant_followup_window_contents(&mut self, ui: &mut egui::Ui) {
        let engine_available = self.engine.is_some();
        let source_seq_id = self.variant_followup_ui.source_seq_id.clone();
        let current_seq_id = self.seq_id.clone().unwrap_or_default();
        let source_missing = !self.variant_followup_sequence_exists(&source_seq_id);
        let has_variant_seed = self.variant_followup_has_variant_seed();
        let promoter_design_mode = if has_variant_seed {
            "variant-seeded"
        } else {
            "gene/promoter-seeded"
        };
        let mut score_track_params_changed = false;
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
                    "Generated promoter windows and TF score tracks are engine-owned evidence; this window only orchestrates the shared operations and displays the portable records.",
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
            ui.small(
                egui::RichText::new(format!("Mode: {promoter_design_mode}"))
                    .color(egui::Color32::from_rgb(100, 116, 139)),
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

            ui.label("TF motifs");
            ui.vertical(|ui| {
                if ui
                    .text_edit_singleline(&mut self.variant_followup_ui.score_track_motifs)
                    .changed()
                {
                    score_track_params_changed = true;
                }
                ui.horizontal_wrapped(|ui| {
                    if ui
                        .small_button("SP1")
                        .on_hover_text("Focus on SP1, useful for the TERT promoter case study.")
                        .clicked()
                    {
                        self.variant_followup_ui.score_track_motifs = "SP1".to_string();
                        score_track_params_changed = true;
                    }
                    if ui
                        .small_button("Yamanaka + NANOG")
                        .on_hover_text("Load POU5F1, SOX2, KLF4, MYC, and NANOG.")
                        .clicked()
                    {
                        self.variant_followup_ui.score_track_motifs =
                            "POU5F1,SOX2,KLF4,MYC,NANOG".to_string();
                        score_track_params_changed = true;
                    }
                    if ui
                        .small_button("TP73 set")
                        .on_hover_text(
                            "Load TP73, SP1, BACH2, and PATZ1 for the TP73 promoter case study.",
                        )
                        .clicked()
                    {
                        self.variant_followup_ui.score_track_motifs =
                            "TP73,SP1,BACH2,PATZ1".to_string();
                        score_track_params_changed = true;
                    }
                });
            });
            ui.end_row();

            ui.label("Score-track range");
            ui.horizontal_wrapped(|ui| {
                if ui
                    .text_edit_singleline(&mut self.variant_followup_ui.score_track_start_0based)
                    .changed()
                {
                    score_track_params_changed = true;
                }
                ui.label("..");
                if ui
                    .text_edit_singleline(
                        &mut self.variant_followup_ui.score_track_end_0based_exclusive,
                    )
                    .changed()
                {
                    score_track_params_changed = true;
                }
                egui::ComboBox::from_id_salt("promoter_design_score_track_value_kind")
                    .selected_text(self.variant_followup_ui.score_track_value_kind.as_str())
                    .show_ui(ui, |ui| {
                        for value_kind in [
                            TfbsScoreTrackValueKind::LlrBits,
                            TfbsScoreTrackValueKind::LlrQuantile,
                            TfbsScoreTrackValueKind::LlrBackgroundQuantile,
                            TfbsScoreTrackValueKind::LlrBackgroundTailLog10,
                            TfbsScoreTrackValueKind::TrueLogOddsBits,
                            TfbsScoreTrackValueKind::TrueLogOddsQuantile,
                            TfbsScoreTrackValueKind::TrueLogOddsBackgroundQuantile,
                            TfbsScoreTrackValueKind::TrueLogOddsBackgroundTailLog10,
                        ] {
                            if ui
                                .selectable_value(
                                    &mut self.variant_followup_ui.score_track_value_kind,
                                    value_kind,
                                    value_kind.as_str(),
                                )
                                .changed()
                            {
                                score_track_params_changed = true;
                            }
                        }
                    });
                if ui
                    .checkbox(
                        &mut self.variant_followup_ui.score_track_clip_negative,
                        "clip negatives",
                    )
                    .changed()
                {
                    score_track_params_changed = true;
                }
            });
            ui.small(
                egui::RichText::new(
                    "Choose raw bits, in-window quantiles, or random-background tail views. Background percentile and -log10 tail plots suppress everything below the 0.95 random-background quantile and are calibrated against a quantized IID random-DNA window model rather than only against a finite sampled background. Negative clipping only affects the raw bit score kinds.",
                )
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
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
        if score_track_params_changed {
            self.variant_followup_ui.cached_score_tracks = None;
        }
        if summary_params_changed || candidate_params_changed {
            self.variant_followup_ui.cached_candidates = None;
        }

        ui.separator();
        ui.horizontal_wrapped(|ui| {
            if ui
                .add_enabled(
                    engine_available && !source_missing && self.tfbs_task.is_none(),
                    egui::Button::new("Show TF score tracks"),
                )
                .on_hover_text(
                    "Compute continuous per-position TF motif scores across the selected promoter range, using positive-only display by default.",
                )
                .clicked()
            {
                self.summarize_variant_followup_score_tracks();
            }
            if ui
                .add_enabled(
                    engine_available && !source_missing,
                    egui::Button::new("Export TF score tracks SVG..."),
                )
                .on_hover_text(
                    "Write the current promoter TF score-track plot as a shared SVG figure through the same engine-owned rendering path used by headless exports.",
                )
                .clicked()
            {
                self.export_variant_followup_score_tracks_svg();
            }
            if ui
                .add_enabled(
                    engine_available && !source_missing,
                    egui::Button::new("Export TF correlation SVG..."),
                )
                .on_hover_text(
                    "Write the current TF synchrony view as a shared SVG heatmap with smoothed/raw Pearson panels and ranked peak-offset pairs.",
                )
                .clicked()
            {
                self.export_variant_followup_score_track_correlation_svg();
            }
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
                    engine_available && !source_missing && has_variant_seed,
                    egui::Button::new("Summarize promoter context"),
                )
                .on_hover_text(
                    "Build the portable variant-promoter context record for ClawBio/OpenClaw handoff when a selected SNP anchors the promoter question.",
                )
                .clicked()
            {
                self.summarize_variant_followup_promoter_context();
            }
            if ui
                .add_enabled(
                    engine_available && !source_missing && has_variant_seed,
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
                .add_enabled(
                    has_variant_seed && has_candidates,
                    egui::Button::new("Extract recommended fragment"),
                )
                .on_hover_text(
                    "Materialize the recommended promoter fragment through the shared ExtractRegion operation.",
                )
                .clicked()
            {
                self.extract_variant_followup_recommended_fragment();
            }
            let has_fragment = !self.variant_followup_ui.fragment_output_id.trim().is_empty();
            if ui
                .add_enabled(
                    has_variant_seed && has_fragment,
                    egui::Button::new("Make reference/alternate inserts"),
                )
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
                .add_enabled(
                    has_variant_seed && has_materialized,
                    egui::Button::new("Preview luciferase pair"),
                )
                .on_hover_text(
                    "Load the pinned local mammalian promoterless luciferase backbone if needed, then build reference/alternate reporter previews.",
                )
                .clicked()
            {
                self.preview_variant_followup_reporter_pair();
            }
            if ui
                .add_enabled(
                    engine_available && !source_missing && has_variant_seed,
                    egui::Button::new("Export handoff bundle"),
                )
                .on_hover_text(
                    "Write a portable ClawBio-facing bundle with promoter-context JSON, reporter-candidate JSON, SVG artifacts, report.md, result.json, and replay commands.",
                )
                .clicked()
            {
                self.export_variant_followup_handoff_bundle();
            }
        });

        ui.separator();
        self.render_variant_followup_score_track_summary(ui);
        ui.add_space(8.0);
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
        let pending_initial_render = self.variant_followup_window_pending_initial_render;
        self.log_promoter_design_status("render begin", pending_initial_render);
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
                                self.render_variant_followup_window_body(
                                    ctx,
                                    ui,
                                    pending_initial_render,
                                );
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
                        self.render_variant_followup_window_body(ctx, ui, pending_initial_render);
                    });
            });

            if crate::app::GENtleApp::viewport_close_requested_or_shortcut(ctx) {
                self.show_variant_followup_window = false;
            }
        });
    }

    fn render_variant_followup_window_body(
        &mut self,
        ctx: &egui::Context,
        ui: &mut egui::Ui,
        pending_initial_render: bool,
    ) {
        if pending_initial_render {
            self.log_promoter_design_status("first frame deferred", true);
            ui.add_space(6.0);
            ui.label(
                egui::RichText::new(
                    "Preparing the Promoter design workspace. Detailed controls will appear on the next repaint.",
                )
                .size(10.0)
                .color(egui::Color32::from_rgb(71, 85, 105)),
            );
            ctx.request_repaint();
            self.variant_followup_window_pending_initial_render = false;
            self.log_promoter_design_status("first frame queued repaint", false);
            return;
        }
        self.log_promoter_design_status("rendering detailed view", false);
        self.render_variant_followup_window_contents(ui);
        self.log_promoter_design_status("render complete", false);
    }
}
