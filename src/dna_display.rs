//! Shared DNA display configuration and visibility policies.

use crate::{
    engine::{
        ConstructReasoningGraph, ConstructRole, EditableStatus, EvidenceClass,
        LinearSequenceLetterLayoutMode, RestrictionEnzymeDisplayMode,
    },
    enzymes::default_preferred_restriction_enzyme_names,
    gc_contents::DEFAULT_SECTION_SIZE_BP,
    restriction_enzyme::{RestrictionEndGeometry, RestrictionEnzymeKey},
};
use std::collections::{BTreeMap, BTreeSet};

use eframe::egui::Color32;

#[derive(Debug, Clone, PartialEq)]
pub struct Selection {
    from: usize,
    to: usize,
    sequence_length: usize,
}

impl Selection {
    pub fn new(from: usize, to: usize, sequence_length: usize) -> Self {
        if sequence_length == 0 {
            return Self {
                from: 0,
                to: 0,
                sequence_length: 0,
            };
        }
        if to > sequence_length {
            Self {
                from: to % sequence_length,
                to: from,
                sequence_length,
            }
        } else {
            Self {
                from,
                to,
                sequence_length,
            }
        }
    }

    pub fn from(&self) -> usize {
        self.from
    }

    pub fn to(&self) -> usize {
        self.to
    }

    pub fn parts(&self) -> Vec<(usize, usize)> {
        if self.sequence_length == 0 {
            return vec![];
        }
        if self.from < self.to {
            vec![(self.from, self.to)]
        } else {
            vec![(self.to, self.sequence_length), (0, self.from)]
        }
    }

    pub fn contains(&self, position: usize) -> bool {
        if self.sequence_length == 0 {
            return false;
        }
        if self.from < self.to {
            position >= self.from && position < self.to
        } else {
            position >= self.from || position < self.to
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum AminoAcidLetters {
    Single,
    Three,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AminoAcidFrame {
    None,
    Automatic,
    Forward(u8),           // 1,2,3
    ReverseCompelment(u8), // 1,2,3
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TfbsDisplayCriteria {
    pub use_llr_bits: bool,
    pub min_llr_bits: f64,
    pub use_llr_quantile: bool,
    pub min_llr_quantile: f64,
    pub use_true_log_odds_bits: bool,
    pub min_true_log_odds_bits: f64,
    pub use_true_log_odds_quantile: bool,
    pub min_true_log_odds_quantile: f64,
}

impl Default for TfbsDisplayCriteria {
    fn default() -> Self {
        Self {
            use_llr_bits: true,
            min_llr_bits: 0.0,
            use_llr_quantile: true,
            min_llr_quantile: 0.95,
            use_true_log_odds_bits: false,
            min_true_log_odds_bits: 0.0,
            use_true_log_odds_quantile: false,
            min_true_log_odds_quantile: 0.95,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct VcfDisplayCriteria {
    pub show_snp: bool,
    pub show_ins: bool,
    pub show_del: bool,
    pub show_sv: bool,
    pub show_other: bool,
    pub pass_only: bool,
    pub use_min_qual: bool,
    pub min_qual: f64,
    pub use_max_qual: bool,
    pub max_qual: f64,
    pub required_info_keys: Vec<String>,
}

impl Default for VcfDisplayCriteria {
    fn default() -> Self {
        Self {
            show_snp: true,
            show_ins: true,
            show_del: true,
            show_sv: true,
            show_other: true,
            pass_only: false,
            use_min_qual: false,
            min_qual: 0.0,
            use_max_qual: false,
            max_qual: 0.0,
            required_info_keys: vec![],
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ConstructReasoningOverlaySpan {
    pub annotation_id: String,
    pub evidence_id: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub strand: Option<String>,
    pub role: ConstructRole,
    pub evidence_class: EvidenceClass,
    pub label: String,
    pub rationale: String,
    pub score: Option<f64>,
    pub confidence: Option<f64>,
    pub context_tags: Vec<String>,
    pub provenance_kind: String,
    pub provenance_refs: Vec<String>,
    pub source_kind: String,
    pub supporting_fact_labels: Vec<String>,
    pub supporting_decision_titles: Vec<String>,
    pub transcript_context_status: Option<String>,
    pub effect_tags: Vec<String>,
    pub editable_status: EditableStatus,
    pub warnings: Vec<String>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ConstructReasoningOverlay {
    pub graph_id: String,
    pub seq_id: String,
    pub objective_title: String,
    pub objective_goal: String,
    pub evidence: Vec<ConstructReasoningOverlaySpan>,
}

impl ConstructReasoningOverlay {
    pub fn from_graph(graph: &ConstructReasoningGraph) -> Self {
        let evidence_by_id = graph
            .evidence
            .iter()
            .map(|row| (row.evidence_id.as_str(), row))
            .collect::<BTreeMap<_, _>>();
        let mut evidence = graph
            .annotation_candidates
            .iter()
            .filter_map(|candidate| {
                let source_evidence = evidence_by_id.get(candidate.evidence_id.as_str())?;
                Some(ConstructReasoningOverlaySpan {
                    annotation_id: candidate.annotation_id.clone(),
                    evidence_id: candidate.evidence_id.clone(),
                    start_0based: candidate.start_0based,
                    end_0based_exclusive: candidate.end_0based_exclusive,
                    strand: candidate.strand.clone(),
                    role: candidate.role,
                    evidence_class: source_evidence.evidence_class,
                    label: if candidate.label.trim().is_empty() {
                        source_evidence.label.clone()
                    } else {
                        candidate.label.clone()
                    },
                    rationale: if candidate.rationale.trim().is_empty() {
                        source_evidence.rationale.clone()
                    } else {
                        candidate.rationale.clone()
                    },
                    score: source_evidence.score,
                    confidence: source_evidence.confidence,
                    context_tags: source_evidence.context_tags.clone(),
                    provenance_kind: source_evidence.provenance_kind.clone(),
                    provenance_refs: source_evidence.provenance_refs.clone(),
                    source_kind: candidate.source_kind.clone(),
                    supporting_fact_labels: candidate.supporting_fact_labels.clone(),
                    supporting_decision_titles: candidate.supporting_decision_titles.clone(),
                    transcript_context_status: candidate.transcript_context_status.clone(),
                    effect_tags: candidate.effect_tags.clone(),
                    editable_status: candidate.editable_status,
                    warnings: candidate.warnings.clone(),
                    notes: candidate.notes.clone(),
                })
            })
            .collect::<Vec<_>>();
        evidence.sort_by(|left, right| {
            left.start_0based
                .cmp(&right.start_0based)
                .then_with(|| left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
                .then_with(|| left.evidence_id.cmp(&right.evidence_id))
        });
        Self {
            graph_id: graph.graph_id.clone(),
            seq_id: graph.seq_id.clone(),
            objective_title: graph.objective.title.clone(),
            objective_goal: graph.objective.goal.clone(),
            evidence,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct UpdateLayoutParts {
    update_map_dna: bool,
    update_map_sequence: bool,
}

impl UpdateLayoutParts {
    pub fn update_all(&mut self) {
        self.update_map_dna = true;
        self.update_map_sequence = true;
    }

    pub fn update_map_dna(&self) -> bool {
        self.update_map_dna
    }

    pub fn update_map_sequence(&self) -> bool {
        self.update_map_sequence
    }

    pub fn map_dna_updated(&mut self) {
        self.update_map_dna = false;
    }

    pub fn map_sequence_updated(&mut self) {
        self.update_map_sequence = false;
    }
}

#[derive(Debug)]
pub struct DnaDisplay {
    show_restriction_enzymes: bool,
    restriction_enzyme_display_mode: RestrictionEnzymeDisplayMode,
    preferred_restriction_enzymes: Vec<String>,
    show_reverse_complement: bool,
    auto_hide_sequence_panel_when_linear_bases_visible: bool,
    sequence_panel_max_text_length_bp: usize,
    show_open_reading_frames: bool,
    suppress_open_reading_frames_for_genome_anchor: bool,
    show_features: bool,
    show_cds_features: bool,
    suppress_cds_features_for_gene_annotations: bool,
    show_gene_features: bool,
    show_mrna_features: bool,
    show_construct_reasoning_overlay: bool,
    construct_reasoning_overlay: Option<ConstructReasoningOverlay>,
    hidden_construct_reasoning_roles: BTreeSet<ConstructRole>,
    hidden_construct_reasoning_evidence_classes: BTreeSet<EvidenceClass>,
    show_contextual_transcript_features: bool,
    show_tfbs: bool,
    regulatory_tracks_near_baseline: bool,
    regulatory_feature_max_view_span_bp: usize,
    hidden_feature_kinds: BTreeSet<String>,
    tfbs_display_criteria: TfbsDisplayCriteria,
    vcf_display_criteria: VcfDisplayCriteria,
    show_gc_contents: bool,
    gc_content_bin_size_bp: usize,
    show_methylation_sites: bool,
    update_layout: UpdateLayoutParts,
    aa_letters: AminoAcidLetters,
    aa_frame: AminoAcidFrame,
    selection: Option<Selection>,
    linear_view_start_bp: usize,
    linear_view_span_bp: usize,
    linear_view_vertical_offset_px: f32,
    linear_sequence_base_text_max_view_span_bp: usize,
    linear_sequence_helical_letters_enabled: bool,
    linear_sequence_helical_max_view_span_bp: usize,
    linear_sequence_condensed_max_view_span_bp: usize,
    linear_sequence_letter_layout_mode: LinearSequenceLetterLayoutMode,
    linear_sequence_helical_phase_offset_bp: usize,
    linear_show_double_strand_bases: bool,
    linear_helical_parallel_strands: bool,
    linear_hide_backbone_when_sequence_bases_visible: bool,
    linear_reverse_strand_use_upside_down_letters: bool,
    reverse_strand_visual_opacity: f32,
    feature_details_font_size: f32,
    linear_external_feature_label_font_size: f32,
    linear_external_feature_label_background_opacity: f32,
}

impl DnaDisplay {
    fn normalized_restriction_enzyme_name(name: &str) -> String {
        name.chars()
            .filter(|c| c.is_ascii_alphanumeric())
            .map(|c| c.to_ascii_uppercase())
            .collect()
    }

    pub fn normalize_preferred_restriction_enzymes(names: &[String]) -> Vec<String> {
        let mut out = Vec::new();
        let mut seen = BTreeSet::new();
        for raw in names {
            let trimmed = raw.trim();
            if trimmed.is_empty() {
                continue;
            }
            let normalized = Self::normalized_restriction_enzyme_name(trimmed);
            if normalized.is_empty() || !seen.insert(normalized) {
                continue;
            }
            out.push(trimmed.to_string());
        }
        out
    }

    pub fn parse_preferred_restriction_enzymes_csv(csv: &str) -> Vec<String> {
        let names = csv
            .split(',')
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>();
        Self::normalize_preferred_restriction_enzymes(&names)
    }

    fn restriction_group_matches_preferred(preferred_names: &[String], names: &[String]) -> bool {
        let preferred = preferred_names
            .iter()
            .map(|name| Self::normalized_restriction_enzyme_name(name))
            .collect::<BTreeSet<_>>();
        names.iter().any(|name| {
            let normalized = Self::normalized_restriction_enzyme_name(name);
            !normalized.is_empty() && preferred.contains(&normalized)
        })
    }

    pub fn restriction_group_matches_mode(
        mode: RestrictionEnzymeDisplayMode,
        preferred_names: &[String],
        key: &RestrictionEnzymeKey,
        names: &[String],
    ) -> bool {
        let is_unique = key.number_of_cuts() == 1;
        let is_preferred = Self::restriction_group_matches_preferred(preferred_names, names);
        match mode {
            RestrictionEnzymeDisplayMode::PreferredOnly => is_preferred,
            RestrictionEnzymeDisplayMode::PreferredAndUnique => is_preferred || is_unique,
            RestrictionEnzymeDisplayMode::UniqueOnly => is_unique,
            RestrictionEnzymeDisplayMode::AllInView => true,
        }
    }

    fn clamp_feature_details_font_size(value: f32) -> f32 {
        value.clamp(8.0, 24.0)
    }

    fn clamp_linear_external_feature_label_font_size(value: f32) -> f32 {
        value.clamp(8.0, 24.0)
    }

    fn clamp_linear_external_feature_label_background_opacity(value: f32) -> f32 {
        value.clamp(0.0, 1.0)
    }

    fn clamp_reverse_strand_visual_opacity(value: f32) -> f32 {
        value.clamp(0.2, 1.0)
    }

    fn clamp_linear_sequence_base_text_max_view_span_bp(value: usize) -> usize {
        value.min(5_000_000)
    }

    fn clamp_sequence_panel_max_text_length_bp(value: usize) -> usize {
        if value == 0 { 0 } else { value.min(5_000_000) }
    }

    fn clamp_gc_content_bin_size_bp(value: usize) -> usize {
        value.clamp(1, 5_000_000)
    }

    fn clamp_linear_sequence_helical_max_view_span_bp(value: usize) -> usize {
        value.min(5_000_000)
    }

    fn clamp_linear_sequence_condensed_max_view_span_bp(value: usize) -> usize {
        value.min(5_000_000)
    }

    fn clamp_linear_view_vertical_offset_px(value: f32) -> f32 {
        if value.is_finite() {
            value.clamp(-10_000.0, 10_000.0)
        } else {
            0.0
        }
    }

    fn clamp_linear_sequence_helical_phase_offset_bp(value: usize) -> usize {
        value % 10
    }

    fn mark_layout_dirty(&mut self) {
        self.update_layout.update_all();
    }

    fn normalize_vcf_info_keys(keys: &[String]) -> Vec<String> {
        let mut dedup = BTreeSet::new();
        for key in keys {
            let normalized = key.trim().to_ascii_uppercase();
            if !normalized.is_empty() {
                dedup.insert(normalized);
            }
        }
        dedup.into_iter().collect()
    }

    pub fn construct_reasoning_overlay_span_visible(
        &self,
        span: &ConstructReasoningOverlaySpan,
    ) -> bool {
        !self.hidden_construct_reasoning_roles.contains(&span.role)
            && !self
                .hidden_construct_reasoning_evidence_classes
                .contains(&span.evidence_class)
    }

    pub fn show_restriction_enzyme_sites(&self) -> bool {
        self.show_restriction_enzymes
    }

    pub fn restriction_enzyme_display_mode(&self) -> RestrictionEnzymeDisplayMode {
        self.restriction_enzyme_display_mode
    }

    pub fn set_restriction_enzyme_display_mode(&mut self, value: RestrictionEnzymeDisplayMode) {
        if self.restriction_enzyme_display_mode != value {
            self.restriction_enzyme_display_mode = value;
            self.mark_layout_dirty();
        }
    }

    pub fn preferred_restriction_enzymes(&self) -> &[String] {
        &self.preferred_restriction_enzymes
    }

    pub fn set_preferred_restriction_enzymes(&mut self, names: Vec<String>) {
        let normalized = Self::normalize_preferred_restriction_enzymes(&names);
        if self.preferred_restriction_enzymes != normalized {
            self.preferred_restriction_enzymes = normalized;
            self.mark_layout_dirty();
        }
    }

    pub fn toggle_show_restriction_enzyme_sites(&mut self) {
        self.show_restriction_enzymes = !self.show_restriction_enzymes;
        self.mark_layout_dirty();
    }

    pub fn set_show_restriction_enzyme_sites(&mut self, value: bool) {
        if self.show_restriction_enzymes != value {
            self.show_restriction_enzymes = value;
            self.mark_layout_dirty();
        }
    }

    pub fn toggle_show_features(&mut self) {
        self.show_features = !self.show_features;
        self.mark_layout_dirty();
    }

    pub fn set_show_features(&mut self, value: bool) {
        if self.show_features != value {
            self.show_features = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_features(&self) -> bool {
        self.show_features
    }

    pub fn show_cds_features(&self) -> bool {
        self.show_cds_features
    }

    pub fn show_cds_features_effective(&self) -> bool {
        self.show_cds_features && !self.suppress_cds_features_for_gene_annotations
    }

    pub fn set_show_cds_features(&mut self, value: bool) {
        if self.show_cds_features != value {
            self.show_cds_features = value;
            self.mark_layout_dirty();
        }
    }

    pub fn suppress_cds_features_for_gene_annotations(&self) -> bool {
        self.suppress_cds_features_for_gene_annotations
    }

    pub fn set_suppress_cds_features_for_gene_annotations(&mut self, value: bool) {
        if self.suppress_cds_features_for_gene_annotations != value {
            self.suppress_cds_features_for_gene_annotations = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_gene_features(&self) -> bool {
        self.show_gene_features
    }

    pub fn set_show_gene_features(&mut self, value: bool) {
        if self.show_gene_features != value {
            self.show_gene_features = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_mrna_features(&self) -> bool {
        self.show_mrna_features
    }

    pub fn set_show_mrna_features(&mut self, value: bool) {
        if self.show_mrna_features != value {
            self.show_mrna_features = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_construct_reasoning_overlay(&self) -> bool {
        self.show_construct_reasoning_overlay
    }

    pub fn toggle_show_construct_reasoning_overlay(&mut self) {
        self.show_construct_reasoning_overlay = !self.show_construct_reasoning_overlay;
        self.mark_layout_dirty();
    }

    pub fn set_show_construct_reasoning_overlay(&mut self, value: bool) {
        if self.show_construct_reasoning_overlay != value {
            self.show_construct_reasoning_overlay = value;
            self.mark_layout_dirty();
        }
    }

    pub fn construct_reasoning_overlay(&self) -> Option<&ConstructReasoningOverlay> {
        self.construct_reasoning_overlay.as_ref()
    }

    pub fn set_construct_reasoning_overlay(&mut self, overlay: Option<ConstructReasoningOverlay>) {
        if self.construct_reasoning_overlay != overlay {
            self.construct_reasoning_overlay = overlay;
            self.mark_layout_dirty();
        }
    }

    pub fn construct_reasoning_role_visible(&self, role: ConstructRole) -> bool {
        !self.hidden_construct_reasoning_roles.contains(&role)
    }

    pub fn hidden_construct_reasoning_roles(&self) -> &BTreeSet<ConstructRole> {
        &self.hidden_construct_reasoning_roles
    }

    pub fn set_construct_reasoning_role_visible(&mut self, role: ConstructRole, visible: bool) {
        let changed = if visible {
            self.hidden_construct_reasoning_roles.remove(&role)
        } else {
            self.hidden_construct_reasoning_roles.insert(role)
        };
        if changed {
            self.mark_layout_dirty();
        }
    }

    pub fn construct_reasoning_evidence_class_visible(&self, class: EvidenceClass) -> bool {
        !self
            .hidden_construct_reasoning_evidence_classes
            .contains(&class)
    }

    pub fn hidden_construct_reasoning_evidence_classes(&self) -> &BTreeSet<EvidenceClass> {
        &self.hidden_construct_reasoning_evidence_classes
    }

    pub fn set_construct_reasoning_evidence_class_visible(
        &mut self,
        class: EvidenceClass,
        visible: bool,
    ) {
        let changed = if visible {
            self.hidden_construct_reasoning_evidence_classes
                .remove(&class)
        } else {
            self.hidden_construct_reasoning_evidence_classes
                .insert(class)
        };
        if changed {
            self.mark_layout_dirty();
        }
    }

    pub fn clear_construct_reasoning_filters(&mut self) {
        if !self.hidden_construct_reasoning_roles.is_empty()
            || !self.hidden_construct_reasoning_evidence_classes.is_empty()
        {
            self.hidden_construct_reasoning_roles.clear();
            self.hidden_construct_reasoning_evidence_classes.clear();
            self.mark_layout_dirty();
        }
    }

    pub fn show_contextual_transcript_features(&self) -> bool {
        self.show_contextual_transcript_features
    }

    pub fn set_show_contextual_transcript_features(&mut self, value: bool) {
        if self.show_contextual_transcript_features != value {
            self.show_contextual_transcript_features = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_tfbs(&self) -> bool {
        self.show_tfbs
    }

    pub fn toggle_show_tfbs(&mut self) {
        self.show_tfbs = !self.show_tfbs;
        self.mark_layout_dirty();
    }

    pub fn set_show_tfbs(&mut self, value: bool) {
        if self.show_tfbs != value {
            self.show_tfbs = value;
            self.mark_layout_dirty();
        }
    }

    pub fn regulatory_tracks_near_baseline(&self) -> bool {
        self.regulatory_tracks_near_baseline
    }

    pub fn set_regulatory_tracks_near_baseline(&mut self, value: bool) {
        if self.regulatory_tracks_near_baseline != value {
            self.regulatory_tracks_near_baseline = value;
            self.mark_layout_dirty();
        }
    }

    pub fn regulatory_feature_max_view_span_bp(&self) -> usize {
        self.regulatory_feature_max_view_span_bp
    }

    pub fn set_regulatory_feature_max_view_span_bp(&mut self, value: usize) {
        if self.regulatory_feature_max_view_span_bp != value {
            self.regulatory_feature_max_view_span_bp = value;
            self.mark_layout_dirty();
        }
    }

    pub fn tfbs_display_criteria(&self) -> TfbsDisplayCriteria {
        self.tfbs_display_criteria
    }

    pub fn set_tfbs_display_criteria(&mut self, criteria: TfbsDisplayCriteria) {
        if self.tfbs_display_criteria != criteria {
            self.tfbs_display_criteria = criteria;
            self.mark_layout_dirty();
        }
    }

    pub fn vcf_display_criteria(&self) -> VcfDisplayCriteria {
        self.vcf_display_criteria.clone()
    }

    pub fn set_vcf_display_criteria(&mut self, mut criteria: VcfDisplayCriteria) {
        criteria.required_info_keys = Self::normalize_vcf_info_keys(&criteria.required_info_keys);
        if self.vcf_display_criteria != criteria {
            self.vcf_display_criteria = criteria;
            self.mark_layout_dirty();
        }
    }

    pub fn hidden_feature_kinds(&self) -> &BTreeSet<String> {
        &self.hidden_feature_kinds
    }

    pub fn set_hidden_feature_kinds(&mut self, mut kinds: BTreeSet<String>) {
        kinds = kinds
            .into_iter()
            .map(|kind| kind.trim().to_ascii_uppercase())
            .filter(|kind| !kind.is_empty())
            .collect();
        if self.hidden_feature_kinds != kinds {
            self.hidden_feature_kinds = kinds;
            self.mark_layout_dirty();
        }
    }

    pub fn feature_kind_visible(&self, kind: &str) -> bool {
        let normalized = kind.trim().to_ascii_uppercase();
        !self.hidden_feature_kinds.contains(&normalized)
    }

    pub fn set_feature_kind_visible(&mut self, kind: &str, visible: bool) {
        let normalized = kind.trim().to_ascii_uppercase();
        if normalized.is_empty() {
            return;
        }
        let changed = if visible {
            self.hidden_feature_kinds.remove(&normalized)
        } else {
            self.hidden_feature_kinds.insert(normalized)
        };
        if changed {
            self.mark_layout_dirty();
        }
    }

    pub fn show_reverse_complement(&self) -> bool {
        self.show_reverse_complement
    }

    pub fn toggle_reverse_complement(&mut self) {
        self.show_reverse_complement = !self.show_reverse_complement;
        self.mark_layout_dirty();
    }

    pub fn set_show_reverse_complement(&mut self, value: bool) {
        if self.show_reverse_complement != value {
            self.show_reverse_complement = value;
            self.mark_layout_dirty();
        }
    }

    pub fn auto_hide_sequence_panel_when_linear_bases_visible(&self) -> bool {
        self.auto_hide_sequence_panel_when_linear_bases_visible
    }

    pub fn set_auto_hide_sequence_panel_when_linear_bases_visible(&mut self, value: bool) {
        if self.auto_hide_sequence_panel_when_linear_bases_visible != value {
            self.auto_hide_sequence_panel_when_linear_bases_visible = value;
            self.mark_layout_dirty();
        }
    }

    pub fn sequence_panel_max_text_length_bp(&self) -> usize {
        Self::clamp_sequence_panel_max_text_length_bp(self.sequence_panel_max_text_length_bp)
    }

    pub fn set_sequence_panel_max_text_length_bp(&mut self, value: usize) {
        let value = Self::clamp_sequence_panel_max_text_length_bp(value);
        if self.sequence_panel_max_text_length_bp != value {
            self.sequence_panel_max_text_length_bp = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_open_reading_frames(&self) -> bool {
        self.show_open_reading_frames
    }

    pub fn show_open_reading_frames_effective(&self) -> bool {
        self.show_open_reading_frames && !self.suppress_open_reading_frames_for_genome_anchor
    }

    pub fn suppress_open_reading_frames_for_genome_anchor(&self) -> bool {
        self.suppress_open_reading_frames_for_genome_anchor
    }

    pub fn set_suppress_open_reading_frames_for_genome_anchor(&mut self, value: bool) {
        if self.suppress_open_reading_frames_for_genome_anchor != value {
            self.suppress_open_reading_frames_for_genome_anchor = value;
            self.mark_layout_dirty();
        }
    }

    pub fn toggle_show_open_reading_frames(&mut self) {
        self.show_open_reading_frames = !self.show_open_reading_frames;
        self.mark_layout_dirty();
    }

    pub fn set_show_open_reading_frames(&mut self, value: bool) {
        if self.show_open_reading_frames != value {
            self.show_open_reading_frames = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_gc_contents(&self) -> bool {
        self.show_gc_contents
    }

    pub fn toggle_show_gc_contents(&mut self) {
        self.show_gc_contents = !self.show_gc_contents;
        self.mark_layout_dirty();
    }

    pub fn set_show_gc_contents(&mut self, value: bool) {
        if self.show_gc_contents != value {
            self.show_gc_contents = value;
            self.mark_layout_dirty();
        }
    }

    pub fn gc_content_bin_size_bp(&self) -> usize {
        Self::clamp_gc_content_bin_size_bp(self.gc_content_bin_size_bp)
    }

    pub fn set_gc_content_bin_size_bp(&mut self, value: usize) {
        let value = Self::clamp_gc_content_bin_size_bp(value);
        if self.gc_content_bin_size_bp != value {
            self.gc_content_bin_size_bp = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_methylation_sites(&self) -> bool {
        self.show_methylation_sites
    }

    pub fn toggle_show_methylation_sites(&mut self) {
        self.show_methylation_sites = !self.show_methylation_sites;
        self.mark_layout_dirty();
    }

    pub fn set_show_methylation_sites(&mut self, value: bool) {
        if self.show_methylation_sites != value {
            self.show_methylation_sites = value;
            self.mark_layout_dirty();
        }
    }

    pub fn update_layout(&self) -> &UpdateLayoutParts {
        &self.update_layout
    }

    pub fn update_layout_mut(&mut self) -> &mut UpdateLayoutParts {
        &mut self.update_layout
    }

    pub fn aa_letters(&self) -> AminoAcidLetters {
        self.aa_letters.clone()
    }

    pub fn set_aa_letters(&mut self, aa_letters: AminoAcidLetters) {
        self.aa_letters = aa_letters;
    }

    pub fn aa_frame(&self) -> AminoAcidFrame {
        self.aa_frame.clone()
    }

    pub fn set_aa_frame(&mut self, aa_frame: AminoAcidFrame) {
        self.aa_frame = aa_frame;
    }

    pub fn selection(&self) -> Option<Selection> {
        self.selection.clone()
    }

    pub fn select(&mut self, selection: Selection) {
        self.selection = Some(selection);
    }

    pub fn deselect(&mut self) {
        self.selection = None;
    }

    pub fn restriction_enzyme_group_color(cuts: usize) -> Color32 {
        match cuts {
            1 => Color32::DARK_RED,
            2 => Color32::DARK_BLUE,
            3 => Color32::DARK_GREEN,
            _ => Color32::BLACK,
        }
    }

    pub fn restriction_enzyme_geometry_color(geometry: RestrictionEndGeometry) -> Color32 {
        match geometry {
            RestrictionEndGeometry::Blunt => Color32::from_rgb(71, 85, 105),
            RestrictionEndGeometry::FivePrimeOverhang(_) => Color32::from_rgb(37, 99, 235),
            RestrictionEndGeometry::ThreePrimeOverhang(_) => Color32::from_rgb(217, 119, 6),
        }
    }

    pub fn linear_view_start_bp(&self) -> usize {
        self.linear_view_start_bp
    }

    pub fn linear_view_span_bp(&self) -> usize {
        self.linear_view_span_bp
    }

    pub fn set_linear_viewport(&mut self, start_bp: usize, span_bp: usize) {
        if self.linear_view_start_bp != start_bp || self.linear_view_span_bp != span_bp {
            self.linear_view_start_bp = start_bp;
            self.linear_view_span_bp = span_bp;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_view_vertical_offset_px(&self) -> f32 {
        Self::clamp_linear_view_vertical_offset_px(self.linear_view_vertical_offset_px)
    }

    pub fn set_linear_view_vertical_offset_px(&mut self, value: f32) {
        let value = Self::clamp_linear_view_vertical_offset_px(value);
        if (self.linear_view_vertical_offset_px - value).abs() > f32::EPSILON {
            self.linear_view_vertical_offset_px = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_sequence_base_text_max_view_span_bp(&self) -> usize {
        Self::clamp_linear_sequence_base_text_max_view_span_bp(
            self.linear_sequence_base_text_max_view_span_bp,
        )
    }

    pub fn set_linear_sequence_base_text_max_view_span_bp(&mut self, value: usize) {
        let value = Self::clamp_linear_sequence_base_text_max_view_span_bp(value);
        if self.linear_sequence_base_text_max_view_span_bp != value {
            self.linear_sequence_base_text_max_view_span_bp = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_sequence_helical_letters_enabled(&self) -> bool {
        self.linear_sequence_helical_letters_enabled
    }

    pub fn set_linear_sequence_helical_letters_enabled(&mut self, value: bool) {
        if self.linear_sequence_helical_letters_enabled != value {
            self.linear_sequence_helical_letters_enabled = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_sequence_helical_max_view_span_bp(&self) -> usize {
        Self::clamp_linear_sequence_helical_max_view_span_bp(
            self.linear_sequence_helical_max_view_span_bp,
        )
    }

    pub fn set_linear_sequence_helical_max_view_span_bp(&mut self, value: usize) {
        let value = Self::clamp_linear_sequence_helical_max_view_span_bp(value);
        if self.linear_sequence_helical_max_view_span_bp != value {
            self.linear_sequence_helical_max_view_span_bp = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_sequence_condensed_max_view_span_bp(&self) -> usize {
        Self::clamp_linear_sequence_condensed_max_view_span_bp(
            self.linear_sequence_condensed_max_view_span_bp,
        )
    }

    pub fn set_linear_sequence_condensed_max_view_span_bp(&mut self, value: usize) {
        let value = Self::clamp_linear_sequence_condensed_max_view_span_bp(value);
        if self.linear_sequence_condensed_max_view_span_bp != value {
            self.linear_sequence_condensed_max_view_span_bp = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_sequence_letter_layout_mode(&self) -> LinearSequenceLetterLayoutMode {
        self.linear_sequence_letter_layout_mode
    }

    pub fn set_linear_sequence_letter_layout_mode(
        &mut self,
        value: LinearSequenceLetterLayoutMode,
    ) {
        if self.linear_sequence_letter_layout_mode != value {
            self.linear_sequence_letter_layout_mode = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_sequence_helical_phase_offset_bp(&self) -> usize {
        Self::clamp_linear_sequence_helical_phase_offset_bp(
            self.linear_sequence_helical_phase_offset_bp,
        )
    }

    pub fn set_linear_sequence_helical_phase_offset_bp(&mut self, value: usize) {
        let value = Self::clamp_linear_sequence_helical_phase_offset_bp(value);
        if self.linear_sequence_helical_phase_offset_bp != value {
            self.linear_sequence_helical_phase_offset_bp = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_show_double_strand_bases(&self) -> bool {
        self.linear_show_double_strand_bases
    }

    pub fn set_linear_show_double_strand_bases(&mut self, value: bool) {
        if self.linear_show_double_strand_bases != value {
            self.linear_show_double_strand_bases = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_helical_parallel_strands(&self) -> bool {
        self.linear_helical_parallel_strands
    }

    pub fn set_linear_helical_parallel_strands(&mut self, value: bool) {
        if self.linear_helical_parallel_strands != value {
            self.linear_helical_parallel_strands = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_hide_backbone_when_sequence_bases_visible(&self) -> bool {
        self.linear_hide_backbone_when_sequence_bases_visible
    }

    pub fn set_linear_hide_backbone_when_sequence_bases_visible(&mut self, value: bool) {
        if self.linear_hide_backbone_when_sequence_bases_visible != value {
            self.linear_hide_backbone_when_sequence_bases_visible = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_reverse_strand_use_upside_down_letters(&self) -> bool {
        self.linear_reverse_strand_use_upside_down_letters
    }

    pub fn set_linear_reverse_strand_use_upside_down_letters(&mut self, value: bool) {
        if self.linear_reverse_strand_use_upside_down_letters != value {
            self.linear_reverse_strand_use_upside_down_letters = value;
            self.mark_layout_dirty();
        }
    }

    pub fn reverse_strand_visual_opacity(&self) -> f32 {
        Self::clamp_reverse_strand_visual_opacity(self.reverse_strand_visual_opacity)
    }

    pub fn set_reverse_strand_visual_opacity(&mut self, value: f32) {
        let value = Self::clamp_reverse_strand_visual_opacity(value);
        if (self.reverse_strand_visual_opacity - value).abs() > f32::EPSILON {
            self.reverse_strand_visual_opacity = value;
            self.mark_layout_dirty();
        }
    }

    pub fn feature_details_font_size(&self) -> f32 {
        Self::clamp_feature_details_font_size(self.feature_details_font_size)
    }

    pub fn set_feature_details_font_size(&mut self, value: f32) {
        let value = Self::clamp_feature_details_font_size(value);
        if (self.feature_details_font_size - value).abs() > f32::EPSILON {
            self.feature_details_font_size = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_external_feature_label_font_size(&self) -> f32 {
        Self::clamp_linear_external_feature_label_font_size(
            self.linear_external_feature_label_font_size,
        )
    }

    pub fn set_linear_external_feature_label_font_size(&mut self, value: f32) {
        let value = Self::clamp_linear_external_feature_label_font_size(value);
        if (self.linear_external_feature_label_font_size - value).abs() > f32::EPSILON {
            self.linear_external_feature_label_font_size = value;
            self.mark_layout_dirty();
        }
    }

    pub fn linear_external_feature_label_background_opacity(&self) -> f32 {
        Self::clamp_linear_external_feature_label_background_opacity(
            self.linear_external_feature_label_background_opacity,
        )
    }

    pub fn set_linear_external_feature_label_background_opacity(&mut self, value: f32) {
        let value = Self::clamp_linear_external_feature_label_background_opacity(value);
        if (self.linear_external_feature_label_background_opacity - value).abs() > f32::EPSILON {
            self.linear_external_feature_label_background_opacity = value;
            self.mark_layout_dirty();
        }
    }
}

impl Default for DnaDisplay {
    fn default() -> Self {
        let mut hidden_feature_kinds = BTreeSet::new();
        hidden_feature_kinds.insert("MISC_FEATURE".to_string());
        Self {
            show_restriction_enzymes: true,
            restriction_enzyme_display_mode: RestrictionEnzymeDisplayMode::default(),
            preferred_restriction_enzymes: default_preferred_restriction_enzyme_names(),
            show_reverse_complement: true,
            auto_hide_sequence_panel_when_linear_bases_visible: false,
            sequence_panel_max_text_length_bp: 200_000,
            show_open_reading_frames: false,
            suppress_open_reading_frames_for_genome_anchor: false,
            show_features: true,
            show_cds_features: true,
            suppress_cds_features_for_gene_annotations: false,
            show_gene_features: true,
            show_mrna_features: true,
            show_construct_reasoning_overlay: true,
            construct_reasoning_overlay: None,
            hidden_construct_reasoning_roles: BTreeSet::new(),
            hidden_construct_reasoning_evidence_classes: BTreeSet::new(),
            show_contextual_transcript_features: true,
            show_tfbs: false,
            regulatory_tracks_near_baseline: false,
            regulatory_feature_max_view_span_bp: 50_000,
            hidden_feature_kinds,
            tfbs_display_criteria: TfbsDisplayCriteria::default(),
            vcf_display_criteria: VcfDisplayCriteria::default(),
            show_gc_contents: true,
            gc_content_bin_size_bp: DEFAULT_SECTION_SIZE_BP,
            show_methylation_sites: false,
            update_layout: UpdateLayoutParts::default(),
            aa_letters: AminoAcidLetters::Single,
            aa_frame: AminoAcidFrame::None,
            selection: None,
            linear_view_start_bp: 0,
            linear_view_span_bp: 0,
            linear_view_vertical_offset_px: 0.0,
            linear_sequence_base_text_max_view_span_bp: 500,
            linear_sequence_helical_letters_enabled: true,
            linear_sequence_helical_max_view_span_bp: 2000,
            linear_sequence_condensed_max_view_span_bp: 1500,
            linear_sequence_letter_layout_mode: LinearSequenceLetterLayoutMode::AutoAdaptive,
            linear_sequence_helical_phase_offset_bp: 0,
            linear_show_double_strand_bases: true,
            linear_helical_parallel_strands: true,
            linear_hide_backbone_when_sequence_bases_visible: false,
            linear_reverse_strand_use_upside_down_letters: true,
            reverse_strand_visual_opacity: 0.55,
            feature_details_font_size: 9.0,
            linear_external_feature_label_font_size: 11.0,
            linear_external_feature_label_background_opacity: 0.9,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{
        AnnotationCandidate, ConstructObjective, DesignDecisionNode, DesignEvidence, DesignFact,
        EvidenceScope,
    };

    #[test]
    fn construct_reasoning_overlay_keeps_generated_and_referenced_sequence_spans() {
        let overlay = ConstructReasoningOverlay::from_graph(&ConstructReasoningGraph {
            graph_id: "graph_demo".to_string(),
            seq_id: "demo".to_string(),
            objective: ConstructObjective {
                title: "Reasoning demo".to_string(),
                goal: "Inspect mapped evidence".to_string(),
                ..ConstructObjective::default()
            },
            evidence: vec![
                DesignEvidence {
                    evidence_id: "host_ctx".to_string(),
                    scope: EvidenceScope::HostProfile,
                    role: ConstructRole::ContextBaggage,
                    evidence_class: EvidenceClass::UserOverride,
                    label: "Propagation host: ecoli_k12".to_string(),
                    rationale: "Objective-selected propagation host".to_string(),
                    ..DesignEvidence::default()
                },
                DesignEvidence {
                    evidence_id: "restriction_site_raw".to_string(),
                    seq_id: "demo".to_string(),
                    scope: EvidenceScope::SequenceSpan,
                    start_0based: 4,
                    end_0based_exclusive: 10,
                    role: ConstructRole::RestrictionSite,
                    evidence_class: EvidenceClass::HardFact,
                    label: "EcoRI".to_string(),
                    rationale: "Raw restriction site evidence".to_string(),
                    ..DesignEvidence::default()
                },
                DesignEvidence {
                    evidence_id: "generated_promoter".to_string(),
                    seq_id: "demo".to_string(),
                    scope: EvidenceScope::SequenceSpan,
                    start_0based: 12,
                    end_0based_exclusive: 48,
                    role: ConstructRole::Promoter,
                    evidence_class: EvidenceClass::ContextEvidence,
                    label: "TP73 promoter window".to_string(),
                    rationale: "Promoter window derived from transcript TSS geometry".to_string(),
                    context_tags: vec!["promoter".to_string(), "generated".to_string()],
                    provenance_kind: "derived_promoter_window".to_string(),
                    ..DesignEvidence::default()
                },
                DesignEvidence {
                    evidence_id: "cds_span".to_string(),
                    seq_id: "demo".to_string(),
                    scope: EvidenceScope::SequenceSpan,
                    start_0based: 55,
                    end_0based_exclusive: 96,
                    role: ConstructRole::Cds,
                    evidence_class: EvidenceClass::ReliableAnnotation,
                    label: "Reporter CDS".to_string(),
                    rationale: "Imported CDS annotation".to_string(),
                    ..DesignEvidence::default()
                },
                DesignEvidence {
                    evidence_id: "variant_span".to_string(),
                    seq_id: "demo".to_string(),
                    scope: EvidenceScope::SequenceSpan,
                    start_0based: 64,
                    end_0based_exclusive: 65,
                    role: ConstructRole::Variant,
                    evidence_class: EvidenceClass::HardFact,
                    label: "rs-demo".to_string(),
                    rationale: "Variant marker".to_string(),
                    ..DesignEvidence::default()
                },
            ],
            facts: vec![DesignFact {
                fact_id: "fact_variant_effect_context".to_string(),
                fact_type: "variant_effect_context".to_string(),
                based_on_evidence_ids: vec!["variant_span".to_string(), "cds_span".to_string()],
                ..DesignFact::default()
            }],
            decisions: vec![DesignDecisionNode {
                decision_id: "decision_variant_effect".to_string(),
                decision_type: "evaluate_variant_effect_context".to_string(),
                input_evidence_ids: vec!["variant_span".to_string(), "cds_span".to_string()],
                ..DesignDecisionNode::default()
            }],
            annotation_candidates: vec![
                AnnotationCandidate {
                    annotation_id: "annotation_generated_promoter".to_string(),
                    evidence_id: "generated_promoter".to_string(),
                    seq_id: "demo".to_string(),
                    start_0based: 12,
                    end_0based_exclusive: 48,
                    strand: Some("+".to_string()),
                    role: ConstructRole::Promoter,
                    label: "TP73 promoter window".to_string(),
                    rationale: "Promoter window derived from transcript TSS geometry".to_string(),
                    source_kind: "generated_annotation".to_string(),
                    ..AnnotationCandidate::default()
                },
                AnnotationCandidate {
                    annotation_id: "annotation_cds_span".to_string(),
                    evidence_id: "cds_span".to_string(),
                    seq_id: "demo".to_string(),
                    start_0based: 55,
                    end_0based_exclusive: 96,
                    strand: Some("+".to_string()),
                    role: ConstructRole::Cds,
                    label: "Reporter CDS".to_string(),
                    rationale: "Imported CDS annotation".to_string(),
                    source_kind: "supporting_annotation".to_string(),
                    supporting_fact_labels: vec!["Variant effect candidates derived".to_string()],
                    supporting_decision_titles: vec!["Evaluate Variant Effect Context".to_string()],
                    ..AnnotationCandidate::default()
                },
                AnnotationCandidate {
                    annotation_id: "annotation_variant_span".to_string(),
                    evidence_id: "variant_span".to_string(),
                    seq_id: "demo".to_string(),
                    start_0based: 64,
                    end_0based_exclusive: 65,
                    strand: Some("+".to_string()),
                    role: ConstructRole::Variant,
                    label: "rs-demo".to_string(),
                    rationale: "Variant marker".to_string(),
                    source_kind: "supporting_annotation".to_string(),
                    supporting_fact_labels: vec!["Variant effect candidates derived".to_string()],
                    supporting_decision_titles: vec!["Evaluate Variant Effect Context".to_string()],
                    transcript_context_status: Some("multi_transcript_ambiguous".to_string()),
                    effect_tags: vec!["promoter_variant_candidate".to_string()],
                    ..AnnotationCandidate::default()
                },
            ],
            ..ConstructReasoningGraph::default()
        });

        assert_eq!(
            overlay
                .evidence
                .iter()
                .map(|row| row.evidence_id.as_str())
                .collect::<Vec<_>>(),
            vec!["generated_promoter", "cds_span", "variant_span"]
        );
    }

    #[test]
    fn construct_reasoning_overlay_keeps_cdna_confirmed_spans_as_annotation_grade() {
        let overlay = ConstructReasoningOverlay::from_graph(&ConstructReasoningGraph {
            graph_id: "graph_demo".to_string(),
            seq_id: "demo".to_string(),
            objective: ConstructObjective {
                title: "Reasoning demo".to_string(),
                goal: "Inspect mapped evidence".to_string(),
                ..ConstructObjective::default()
            },
            evidence: vec![DesignEvidence {
                evidence_id: "confirmed_exon".to_string(),
                seq_id: "demo".to_string(),
                scope: EvidenceScope::SequenceSpan,
                start_0based: 120,
                end_0based_exclusive: 180,
                role: ConstructRole::Exon,
                evidence_class: EvidenceClass::HardFact,
                label: "Confirmed exon".to_string(),
                rationale: "cDNA-confirmed exon annotation".to_string(),
                context_tags: vec!["exon".to_string(), "cdna_confirmed".to_string()],
                ..DesignEvidence::default()
            }],
            annotation_candidates: vec![AnnotationCandidate {
                annotation_id: "annotation_confirmed_exon".to_string(),
                evidence_id: "confirmed_exon".to_string(),
                seq_id: "demo".to_string(),
                start_0based: 120,
                end_0based_exclusive: 180,
                strand: Some("+".to_string()),
                role: ConstructRole::Exon,
                label: "Confirmed exon".to_string(),
                rationale: "cDNA-confirmed exon annotation".to_string(),
                source_kind: "confirmed_annotation".to_string(),
                ..AnnotationCandidate::default()
            }],
            ..ConstructReasoningGraph::default()
        });

        assert_eq!(overlay.evidence.len(), 1);
        assert_eq!(overlay.evidence[0].evidence_id, "confirmed_exon");
    }

    #[test]
    fn construct_reasoning_overlay_keeps_written_back_generated_promoter_features() {
        let overlay = ConstructReasoningOverlay::from_graph(&ConstructReasoningGraph {
            graph_id: "graph_demo".to_string(),
            seq_id: "demo".to_string(),
            objective: ConstructObjective {
                title: "Reasoning demo".to_string(),
                goal: "Inspect mapped evidence".to_string(),
                ..ConstructObjective::default()
            },
            evidence: vec![DesignEvidence {
                evidence_id: "annotated_promoter".to_string(),
                seq_id: "demo".to_string(),
                scope: EvidenceScope::SequenceSpan,
                start_0based: 40,
                end_0based_exclusive: 120,
                role: ConstructRole::Promoter,
                evidence_class: EvidenceClass::ReliableAnnotation,
                label: "Generated promoter feature".to_string(),
                rationale: "Promoter feature written back by AnnotatePromoterWindows".to_string(),
                context_tags: vec![
                    "promoter".to_string(),
                    "annotate_promoter_windows".to_string(),
                ],
                provenance_kind: "sequence_feature_annotation".to_string(),
                ..DesignEvidence::default()
            }],
            annotation_candidates: vec![AnnotationCandidate {
                annotation_id: "annotation_written_back_promoter".to_string(),
                evidence_id: "annotated_promoter".to_string(),
                seq_id: "demo".to_string(),
                start_0based: 40,
                end_0based_exclusive: 120,
                strand: Some("+".to_string()),
                role: ConstructRole::Promoter,
                label: "Generated promoter feature".to_string(),
                rationale: "Promoter feature written back by AnnotatePromoterWindows".to_string(),
                source_kind: "generated_annotation".to_string(),
                ..AnnotationCandidate::default()
            }],
            ..ConstructReasoningGraph::default()
        });

        assert_eq!(overlay.evidence.len(), 1);
        assert_eq!(overlay.evidence[0].evidence_id, "annotated_promoter");
    }
}
