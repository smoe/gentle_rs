use crate::gc_contents::DEFAULT_SECTION_SIZE_BP;
use std::collections::BTreeSet;

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
    show_reverse_complement: bool,
    auto_hide_sequence_panel_when_linear_bases_visible: bool,
    show_open_reading_frames: bool,
    suppress_open_reading_frames_for_genome_anchor: bool,
    show_features: bool,
    show_cds_features: bool,
    suppress_cds_features_for_gene_annotations: bool,
    show_gene_features: bool,
    show_mrna_features: bool,
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
    linear_sequence_base_text_max_view_span_bp: usize,
    linear_sequence_helical_letters_enabled: bool,
    linear_sequence_helical_max_view_span_bp: usize,
    linear_show_double_strand_bases: bool,
    linear_hide_backbone_when_sequence_bases_visible: bool,
    linear_reverse_strand_use_upside_down_letters: bool,
    feature_details_font_size: f32,
    linear_external_feature_label_font_size: f32,
    linear_external_feature_label_background_opacity: f32,
}

impl DnaDisplay {
    fn clamp_feature_details_font_size(value: f32) -> f32 {
        value.clamp(8.0, 24.0)
    }

    fn clamp_linear_external_feature_label_font_size(value: f32) -> f32 {
        value.clamp(8.0, 24.0)
    }

    fn clamp_linear_external_feature_label_background_opacity(value: f32) -> f32 {
        value.clamp(0.0, 1.0)
    }

    fn clamp_linear_sequence_base_text_max_view_span_bp(value: usize) -> usize {
        value.min(5_000_000)
    }

    fn clamp_gc_content_bin_size_bp(value: usize) -> usize {
        value.clamp(1, 5_000_000)
    }

    fn clamp_linear_sequence_helical_max_view_span_bp(value: usize) -> usize {
        value.min(5_000_000)
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

    pub fn show_restriction_enzyme_sites(&self) -> bool {
        self.show_restriction_enzymes
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

    pub fn linear_show_double_strand_bases(&self) -> bool {
        self.linear_show_double_strand_bases
    }

    pub fn set_linear_show_double_strand_bases(&mut self, value: bool) {
        if self.linear_show_double_strand_bases != value {
            self.linear_show_double_strand_bases = value;
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
            show_reverse_complement: true,
            auto_hide_sequence_panel_when_linear_bases_visible: false,
            show_open_reading_frames: false,
            suppress_open_reading_frames_for_genome_anchor: false,
            show_features: true,
            show_cds_features: true,
            suppress_cds_features_for_gene_annotations: false,
            show_gene_features: true,
            show_mrna_features: true,
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
            linear_sequence_base_text_max_view_span_bp: 500,
            linear_sequence_helical_letters_enabled: false,
            linear_sequence_helical_max_view_span_bp: 2000,
            linear_show_double_strand_bases: true,
            linear_hide_backbone_when_sequence_bases_visible: false,
            linear_reverse_strand_use_upside_down_letters: true,
            feature_details_font_size: 8.25,
            linear_external_feature_label_font_size: 11.0,
            linear_external_feature_label_background_opacity: 0.9,
        }
    }
}
