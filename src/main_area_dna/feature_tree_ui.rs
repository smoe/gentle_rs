//! Feature-tree model, filtering, grouping, and side-panel rendering for `MainAreaDna`.
//!
//! This submodule keeps the sequence-window feature browser close to its cache
//! and grouping helpers while preserving the parent `MainAreaDna` orchestration
//! state.

use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub(super) enum FeatureTreeGroupingMode {
    #[default]
    Off,
    Auto,
    Always,
}

impl FeatureTreeGroupingMode {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Off => "Off",
            Self::Auto => "Auto (duplicates)",
            Self::Always => "Always",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum FeatureTreeFocusPreset {
    Cloning,
    Regulatory,
    Repeats,
    Tracks,
}

impl FeatureTreeFocusPreset {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Cloning => "Cloning",
            Self::Regulatory => "Regulatory",
            Self::Repeats => "Repeats",
            Self::Tracks => "Tracks",
        }
    }

    pub(super) fn filter_text(self) -> &'static str {
        match self {
            Self::Cloning => "",
            Self::Regulatory => "regulatory",
            Self::Repeats => "repeat",
            Self::Tracks => "track",
        }
    }

    pub(super) fn grouping_mode(self) -> FeatureTreeGroupingMode {
        match self {
            Self::Cloning => FeatureTreeGroupingMode::Auto,
            Self::Regulatory | Self::Repeats | Self::Tracks => FeatureTreeGroupingMode::Always,
        }
    }

    pub(super) fn hover_text(self) -> &'static str {
        match self {
            Self::Cloning => {
                "Clear the tree filter and use automatic grouping for a cloning-oriented feature list"
            }
            Self::Regulatory => {
                "Focus the tree on regulatory annotations and keep nested grouping openable"
            }
            Self::Repeats => {
                "Focus the tree on repeat/rmsk annotations grouped by repeat class and family"
            }
            Self::Tracks => {
                "Focus the tree on imported track and array annotations grouped by source/experiment"
            }
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct RegulatoryFeatureGrouping {
    pub(super) primary_key: String,
    pub(super) primary_label: String,
    pub(super) secondary_key: Option<String>,
    pub(super) secondary_label: Option<String>,
    pub(super) display_label: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct FeatureTreeEntry {
    pub(super) id: usize,
    pub(super) feature_label_full: String,
    pub(super) feature_label: String,
    pub(super) range_label: String,
    pub(super) subgroup_key: Option<String>,
    pub(super) subgroup_label: Option<String>,
    pub(super) prefer_grouped_label: bool,
    pub(super) show_range_inline_when_ungrouped: bool,
    pub(super) visible_in_view: bool,
    pub(super) is_regulatory: bool,
    pub(super) is_repeat: bool,
    pub(super) is_track: bool,
    pub(super) is_array_track: bool,
    pub(super) is_tfbs: bool,
    pub(super) disable_grouping: bool,
    pub(super) supports_splicing_expert: bool,
    pub(super) supports_variant_followup: bool,
    pub(super) can_seed_promoter_anchor: bool,
    pub(super) regulatory_primary_group_key: Option<String>,
    pub(super) regulatory_primary_group_label: Option<String>,
    pub(super) regulatory_secondary_group_key: Option<String>,
    pub(super) regulatory_secondary_group_label: Option<String>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct FeatureTreeLayerSummaryBucket {
    pub(super) visible_count: usize,
    pub(super) total_count: usize,
}

impl FeatureTreeLayerSummaryBucket {
    pub(super) fn add(&mut self, visible_in_view: bool) {
        self.total_count = self.total_count.saturating_add(1);
        if visible_in_view {
            self.visible_count = self.visible_count.saturating_add(1);
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct FeatureTreeLayerSummary {
    pub(super) core: FeatureTreeLayerSummaryBucket,
    pub(super) regulatory: FeatureTreeLayerSummaryBucket,
    pub(super) repeats: FeatureTreeLayerSummaryBucket,
    pub(super) tracks: FeatureTreeLayerSummaryBucket,
    pub(super) arrays: FeatureTreeLayerSummaryBucket,
    pub(super) tfbs: FeatureTreeLayerSummaryBucket,
    pub(super) other: FeatureTreeLayerSummaryBucket,
}

impl FeatureTreeLayerSummary {
    pub(super) fn add_entry(&mut self, kind: &str, entry: &FeatureTreeEntry) {
        let bucket = if entry.is_array_track {
            &mut self.arrays
        } else if entry.is_track {
            &mut self.tracks
        } else if entry.is_repeat {
            &mut self.repeats
        } else if entry.is_tfbs {
            &mut self.tfbs
        } else if entry.is_regulatory {
            &mut self.regulatory
        } else {
            match kind.trim().to_ascii_uppercase().as_str() {
                "CDS" | "GENE" | "MRNA" => &mut self.core,
                _ => &mut self.other,
            }
        };
        bucket.add(entry.visible_in_view);
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum FeatureCopyPayloadKind {
    Identifier,
    Description,
    PopupText,
}

impl FeatureCopyPayloadKind {
    pub(super) fn menu_label(self) -> &'static str {
        match self {
            Self::Identifier => "Copy identifier",
            Self::Description => "Copy description",
            Self::PopupText => "Copy tooltip text",
        }
    }

    pub(super) fn hover_text(self) -> &'static str {
        match self {
            Self::Identifier => {
                "Copy the most specific stable feature identifier, such as transcript_id or protein_id"
            }
            Self::Description => {
                "Copy the human-readable feature description/title shown by GENtle"
            }
            Self::PopupText => {
                "Copy the full feature popup text: title, range, and qualifier details"
            }
        }
    }

    pub(super) fn status_label(self) -> &'static str {
        match self {
            Self::Identifier => "identifier",
            Self::Description => "description",
            Self::PopupText => "tooltip text",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct FeatureTreeNamedGroup {
    pub(super) key: String,
    pub(super) label: String,
    pub(super) entry_indices: Vec<usize>,
    pub(super) visible_count: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct FeatureTreeRegulatoryPrimaryGroup {
    pub(super) key: String,
    pub(super) label: String,
    pub(super) entry_indices: Vec<usize>,
    pub(super) visible_count: usize,
    pub(super) secondary_groups: Vec<FeatureTreeNamedGroup>,
    pub(super) ungrouped_entry_indices: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct FeatureTreeKindGroup {
    pub(super) kind: String,
    pub(super) entries: Vec<FeatureTreeEntry>,
    pub(super) visible_count: usize,
    pub(super) grouped_entries: Vec<FeatureTreeNamedGroup>,
    pub(super) regulatory_primary_groups: Vec<FeatureTreeRegulatoryPrimaryGroup>,
    pub(super) ungrouped_entry_indices: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct FeatureTreeCacheKey {
    pub(super) seq_len: usize,
    pub(super) feature_count: usize,
    pub(super) viewport: Option<(usize, usize)>,
    pub(super) grouping_mode: FeatureTreeGroupingMode,
    pub(super) feature_filter_text: String,
    pub(super) show_cds_features: bool,
    pub(super) show_gene_features: bool,
    pub(super) show_mrna_features: bool,
    pub(super) show_repeat_features: bool,
    pub(super) show_array_features: bool,
    pub(super) show_contextual_transcript_features: bool,
    pub(super) show_tfbs: bool,
    pub(super) tfbs_display_criteria: TfbsDisplayCriteria,
    pub(super) vcf_display_criteria: VcfDisplayCriteria,
    pub(super) hidden_feature_kinds: BTreeSet<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct FeatureTreeComputedModel {
    pub(super) filter_total_count: usize,
    pub(super) filter_matched_count: usize,
    pub(super) groups: Vec<FeatureTreeKindGroup>,
}

#[derive(Clone, Debug)]
pub(super) struct FeatureTreeCache {
    pub(super) key: FeatureTreeCacheKey,
    pub(super) model: FeatureTreeComputedModel,
}

impl MainAreaDna {
    pub(super) fn format_feature_tree_count_label(
        label: &str,
        visible_count: usize,
        total_count: usize,
        viewport_limited: bool,
    ) -> String {
        if viewport_limited && visible_count != total_count {
            format!("{label} ({visible_count}/{total_count})")
        } else {
            format!("{label} ({total_count})")
        }
    }

    pub(super) fn format_feature_tree_layer_summary_label(
        label: &str,
        bucket: FeatureTreeLayerSummaryBucket,
        viewport_limited: bool,
    ) -> Option<String> {
        if bucket.total_count == 0 {
            return None;
        }
        Some(Self::format_feature_tree_count_label(
            label,
            bucket.visible_count,
            bucket.total_count,
            viewport_limited,
        ))
    }

    pub(super) fn feature_tree_layer_summary(
        model: &FeatureTreeComputedModel,
    ) -> FeatureTreeLayerSummary {
        let mut summary = FeatureTreeLayerSummary::default();
        for group in &model.groups {
            for entry in &group.entries {
                summary.add_entry(&group.kind, entry);
            }
        }
        summary
    }

    pub(super) fn feature_tree_layer_summary_labels(
        model: &FeatureTreeComputedModel,
        viewport_limited: bool,
    ) -> Vec<String> {
        let summary = Self::feature_tree_layer_summary(model);
        [
            ("core", summary.core),
            ("regulatory", summary.regulatory),
            ("repeats", summary.repeats),
            ("tracks", summary.tracks),
            ("array", summary.arrays),
            ("TFBS", summary.tfbs),
            ("other", summary.other),
        ]
        .into_iter()
        .filter_map(|(label, bucket)| {
            Self::format_feature_tree_layer_summary_label(label, bucket, viewport_limited)
        })
        .collect()
    }

    pub(super) fn apply_feature_tree_focus_preset(&mut self, preset: FeatureTreeFocusPreset) {
        self.feature_tree_filter = preset.filter_text().to_string();
        self.feature_tree_grouping_mode = preset.grouping_mode();
        self.pending_feature_tree_scroll_to = None;
        self.save_engine_ops_state();
    }

    pub(super) fn feature_tree_subgroup_label(
        feature: &gb_io::seq::Feature,
        feature_label: &str,
        grouping_mode: FeatureTreeGroupingMode,
    ) -> Option<String> {
        if matches!(grouping_mode, FeatureTreeGroupingMode::Off) {
            return None;
        }
        if RenderDna::is_tfbs_feature(feature) {
            return RenderDna::tfbs_group_label(feature);
        }
        if RenderDna::is_restriction_site_feature(feature) {
            return RenderDna::restriction_site_group_label(feature);
        }
        if RenderDna::is_repeat_feature(feature) {
            return RenderDna::repeat_group_label(feature);
        }
        if RenderDna::is_track_feature(feature) {
            return RenderDna::track_group_label(feature);
        }
        if feature.kind.to_string().eq_ignore_ascii_case("GENE") {
            return None;
        }
        if RenderDna::is_regulatory_feature(feature)
            && !RenderDna::is_track_feature(feature)
            && !RenderDna::is_tfbs_feature(feature)
        {
            return None;
        }
        if feature.kind.to_string().eq_ignore_ascii_case("MRNA")
            && let Some(gene_label) = Self::feature_tree_first_nonempty_qualifier(
                feature,
                &["gene", "gene_name", "locus_tag", "gene_id"],
            )
        {
            return Some(gene_label);
        }
        let normalized_label = feature_label.trim();
        if normalized_label.is_empty() {
            None
        } else {
            Some(normalized_label.to_string())
        }
    }

    pub(super) fn feature_tree_should_group(
        grouping_mode: FeatureTreeGroupingMode,
        subgroup_count: usize,
    ) -> bool {
        match grouping_mode {
            FeatureTreeGroupingMode::Off => false,
            FeatureTreeGroupingMode::Auto => subgroup_count > 1,
            FeatureTreeGroupingMode::Always => true,
        }
    }

    pub(super) fn feature_tree_first_nonempty_qualifier(
        feature: &gb_io::seq::Feature,
        keys: &[&str],
    ) -> Option<String> {
        for key in keys {
            for value in feature.qualifier_values(key) {
                let normalized = value.split_whitespace().collect::<Vec<_>>().join(" ");
                let normalized = normalized.trim().to_string();
                if !normalized.is_empty() {
                    return Some(normalized);
                }
            }
        }
        None
    }

    pub(super) fn feature_copy_identifier_text(feature: &gb_io::seq::Feature) -> String {
        Self::feature_tree_first_nonempty_qualifier(
            feature,
            &[
                "transcript_id",
                "mrna_id",
                "rna_id",
                "protein_id",
                "gene",
                "gene_id",
                "locus_tag",
                "standard_name",
                "label",
                "name",
                "db_xref",
            ],
        )
        .unwrap_or_else(|| RenderDna::feature_name(feature))
    }

    pub(super) fn feature_copy_description_text(feature: &gb_io::seq::Feature) -> String {
        Self::feature_tree_first_nonempty_qualifier(
            feature,
            &[
                "product",
                "label",
                "name",
                "standard_name",
                "note",
                "function",
                "gene",
                "gene_synonym",
                "bound_moiety",
                "regulatory_class",
            ],
        )
        .unwrap_or_else(|| RenderDna::feature_name(feature))
    }

    pub(super) fn feature_copy_popup_text(feature: &gb_io::seq::Feature) -> String {
        let mut lines = Vec::new();
        let title = RenderDna::feature_name(feature);
        if !title.trim().is_empty() {
            lines.push(title);
        }
        let range = RenderDna::feature_range_text(feature);
        if !range.trim().is_empty() {
            lines.push(range);
        }
        lines.extend(RenderDna::feature_detail_lines(feature));
        if lines.is_empty() {
            lines.push(feature.kind.to_string());
        }
        lines.join("\n")
    }

    pub(super) fn feature_copy_payload(
        feature: &gb_io::seq::Feature,
        kind: FeatureCopyPayloadKind,
    ) -> String {
        match kind {
            FeatureCopyPayloadKind::Identifier => Self::feature_copy_identifier_text(feature),
            FeatureCopyPayloadKind::Description => Self::feature_copy_description_text(feature),
            FeatureCopyPayloadKind::PopupText => Self::feature_copy_popup_text(feature),
        }
    }

    pub(super) fn copy_feature_payload_to_clipboard(
        &mut self,
        ctx: &egui::Context,
        feature_id: usize,
        kind: FeatureCopyPayloadKind,
    ) {
        let payload = {
            let dna = self.dna.read().expect("DNA lock poisoned");
            dna.features()
                .get(feature_id)
                .map(|feature| Self::feature_copy_payload(feature, kind))
        };
        match payload {
            Some(payload) if !payload.trim().is_empty() => {
                ctx.copy_text(payload);
                self.op_status = format!(
                    "Copied feature {} for feature #{}",
                    kind.status_label(),
                    feature_id + 1
                );
            }
            Some(_) => {
                self.op_status = format!(
                    "Feature #{} has no copyable {}",
                    feature_id + 1,
                    kind.status_label()
                );
            }
            None => {
                self.op_status = format!("Feature #{} is no longer available", feature_id + 1);
            }
        }
    }

    pub(super) fn render_feature_copy_context_menu_items(
        &mut self,
        ui: &mut egui::Ui,
        feature_id: usize,
    ) -> bool {
        let mut copied = false;
        for kind in [
            FeatureCopyPayloadKind::Identifier,
            FeatureCopyPayloadKind::Description,
            FeatureCopyPayloadKind::PopupText,
        ] {
            if ui
                .button(kind.menu_label())
                .on_hover_text(kind.hover_text())
                .clicked()
            {
                self.copy_feature_payload_to_clipboard(ui.ctx(), feature_id, kind);
                copied = true;
            }
        }
        if copied {
            ui.close();
        }
        copied
    }

    pub(super) fn feature_tree_rna_discriminator(
        feature: &gb_io::seq::Feature,
        range_label: &str,
    ) -> Option<String> {
        let kind = feature.kind.to_string().to_ascii_uppercase();
        if !kind.contains("RNA") {
            return None;
        }
        Self::feature_tree_first_nonempty_qualifier(
            feature,
            &["transcript_id", "mrna_id", "rna_id"],
        )
        .or_else(|| {
            Self::feature_tree_first_nonempty_qualifier(feature, &["transcript_variant", "variant"])
        })
        .or_else(|| {
            let trimmed = range_label.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        })
    }

    pub(super) fn feature_tree_display_label(
        feature: &gb_io::seq::Feature,
        base_label: &str,
        range_label: &str,
    ) -> String {
        let base_label = base_label.trim();
        if base_label.is_empty() {
            return String::new();
        }
        let Some(discriminator) = Self::feature_tree_rna_discriminator(feature, range_label) else {
            return base_label.to_string();
        };
        let discriminator = discriminator.trim();
        if discriminator.is_empty() {
            return base_label.to_string();
        }
        if base_label
            .to_ascii_lowercase()
            .contains(&discriminator.to_ascii_lowercase())
        {
            return base_label.to_string();
        }
        format!("{base_label} [{discriminator}]")
    }

    pub(super) fn trim_feature_tree_prefix_separators(value: &str) -> String {
        value
            .trim_start_matches(|ch: char| {
                ch.is_whitespace() || matches!(ch, ':' | '-' | '_' | ';' | ',')
            })
            .trim()
            .to_string()
    }

    pub(super) fn strip_case_insensitive_prefix(value: &str, prefix: &str) -> Option<String> {
        let trimmed = value.trim_start();
        let lower_trimmed = trimmed.to_ascii_lowercase();
        let lower_prefix = prefix.trim().to_ascii_lowercase();
        if lower_prefix.is_empty() {
            return None;
        }
        if lower_trimmed == lower_prefix {
            return Some(String::new());
        }
        for separator in [":", " ", "_", "-"] {
            let pattern = format!("{lower_prefix}{separator}");
            if lower_trimmed.starts_with(&pattern) {
                return Some(trimmed[pattern.len()..].to_string());
            }
        }
        None
    }

    pub(super) fn extract_histone_like_marker(value: &str) -> Option<String> {
        let token = value
            .split_whitespace()
            .next()
            .unwrap_or_default()
            .trim_matches(|ch: char| matches!(ch, ':' | ';' | ',' | '(' | ')' | '[' | ']'));
        if token.len() < 4 {
            return None;
        }
        let starts_with_h = token
            .chars()
            .next()
            .map(|ch| ch.eq_ignore_ascii_case(&'h'))
            .unwrap_or(false);
        if !starts_with_h {
            return None;
        }
        let mut has_digit = false;
        for ch in token.chars() {
            if ch.is_ascii_digit() {
                has_digit = true;
            }
            if !(ch.is_ascii_alphanumeric() || ch == '-' || ch == '_') {
                return None;
            }
        }
        if has_digit {
            Some(token.to_string())
        } else {
            None
        }
    }

    pub(super) fn is_regulatory_marker_stopword(word: &str) -> bool {
        matches!(
            word,
            "active"
                | "region"
                | "enhancer"
                | "silencer"
                | "promoter"
                | "fragment"
                | "regulatory"
                | "putative"
                | "candidate"
                | "sequence"
                | "assembly"
                | "coordinates"
                | "cell"
                | "cells"
                | "line"
                | "human"
                | "mouse"
                | "hesc"
                | "atac"
                | "chip"
                | "starr"
                | "seq"
                | "grch37"
                | "grch38"
                | "hg19"
                | "hg38"
        )
    }

    pub(super) fn looks_like_protein_marker_token(token: &str) -> bool {
        if token.len() < 2 {
            return false;
        }
        let mut has_alpha = false;
        let mut has_digit = false;
        let mut letters = Vec::new();
        for ch in token.chars() {
            if ch.is_ascii_alphabetic() {
                has_alpha = true;
                letters.push(ch);
            } else if ch.is_ascii_digit() {
                has_digit = true;
            } else if ch == '_' {
                continue;
            } else {
                return false;
            }
        }
        if !has_alpha {
            return false;
        }
        let uppercase_letters = letters.iter().filter(|ch| ch.is_ascii_uppercase()).count();
        let lowercase_letters = letters.iter().filter(|ch| ch.is_ascii_lowercase()).count();
        let uppercase_like = uppercase_letters >= 2 && lowercase_letters == 0;
        has_digit || uppercase_like
    }

    pub(super) fn extract_regulatory_marker_token(value: &str) -> Option<String> {
        let token = value
            .split_whitespace()
            .next()
            .unwrap_or_default()
            .trim_matches(|ch: char| matches!(ch, ':' | ';' | ',' | '(' | ')' | '[' | ']'));
        if token.is_empty() {
            return None;
        }
        if let Some(histone) = Self::extract_histone_like_marker(token) {
            return Some(histone);
        }

        let mut candidate = token.to_string();
        for separator in ['-', '/', '|'] {
            if let Some((head, _)) = candidate.split_once(separator) {
                let head = head.trim();
                if !head.is_empty() {
                    candidate = head.to_string();
                    break;
                }
            }
        }
        let candidate = candidate
            .trim_matches(|ch: char| matches!(ch, ':' | ';' | ',' | '(' | ')' | '[' | ']'));
        if candidate.is_empty() {
            return None;
        }

        let lower = candidate.to_ascii_lowercase();
        if Self::is_regulatory_marker_stopword(&lower) {
            return None;
        }
        if lower.starts_with("chr")
            && lower
                .chars()
                .nth(3)
                .map(|ch| ch.is_ascii_digit())
                .unwrap_or(false)
        {
            return None;
        }
        if Self::looks_like_protein_marker_token(candidate) {
            Some(candidate.to_string())
        } else {
            None
        }
    }

    pub(super) fn derive_regulatory_feature_grouping(
        feature: &gb_io::seq::Feature,
        full_label: &str,
    ) -> Option<RegulatoryFeatureGrouping> {
        if !RenderDna::is_regulatory_feature(feature)
            || RenderDna::is_track_feature(feature)
            || RenderDna::is_tfbs_feature(feature)
        {
            return None;
        }

        let mut detected_primary =
            Self::feature_tree_first_nonempty_qualifier(feature, &["regulatory_class"])
                .map(|value| value.to_ascii_lowercase());
        if detected_primary.is_none() {
            let label_lower = full_label.trim().to_ascii_lowercase();
            for candidate in ["enhancer", "silencer"] {
                if label_lower == candidate
                    || label_lower.starts_with(&format!("{candidate}:"))
                    || label_lower.starts_with(&format!("{candidate} "))
                    || label_lower.starts_with(&format!("{candidate}_"))
                    || label_lower.starts_with(&format!("{candidate}-"))
                {
                    detected_primary = Some(candidate.to_string());
                    break;
                }
            }
        }

        let primary = match detected_primary.as_deref() {
            Some("enhancer") => "enhancer".to_string(),
            Some("silencer") => "silencer".to_string(),
            _ => "other".to_string(),
        };

        let mut display_label = if primary == "other" {
            full_label.to_string()
        } else {
            Self::strip_case_insensitive_prefix(full_label, &primary)
                .unwrap_or_else(|| full_label.to_string())
        };
        display_label = Self::trim_feature_tree_prefix_separators(&display_label);

        let mut secondary_key: Option<String> = None;
        let mut secondary_label: Option<String> = None;

        if let Some(stripped_active) =
            Self::strip_case_insensitive_prefix(&display_label, "active region")
                .or_else(|| Self::strip_case_insensitive_prefix(&display_label, "active_region"))
        {
            secondary_key = Some("active_region".to_string());
            secondary_label = Some("active region".to_string());
            let candidate = Self::trim_feature_tree_prefix_separators(&stripped_active);
            if candidate.chars().count() >= 2 {
                display_label = candidate;
            }
        } else if primary == "enhancer"
            && let Some(marker) = Self::extract_regulatory_marker_token(&display_label)
        {
            secondary_key = Some(format!("marker:{}", marker.to_ascii_lowercase()));
            secondary_label = Some(marker.clone());
            if let Some(after_marker) = Self::strip_case_insensitive_prefix(&display_label, &marker)
            {
                let candidate = Self::trim_feature_tree_prefix_separators(&after_marker);
                if candidate.chars().count() >= 2 {
                    display_label = candidate;
                }
            }
        }

        if display_label.trim().is_empty() {
            display_label = full_label.trim().to_string();
        }

        Some(RegulatoryFeatureGrouping {
            primary_key: primary.clone(),
            primary_label: primary,
            secondary_key,
            secondary_label,
            display_label,
        })
    }

    pub(super) fn feature_tree_matches_filter(
        feature: &gb_io::seq::Feature,
        filter_text: &str,
        kind_label: &str,
        display_label: &str,
        range_label: &str,
    ) -> bool {
        let terms = Self::parse_feature_tree_filter_terms(filter_text);
        if terms.is_empty() {
            return true;
        }

        let all_terms = {
            let mut values = vec![
                kind_label.to_string(),
                display_label.to_string(),
                range_label.to_string(),
                feature.kind.to_string(),
                RenderDna::feature_name(feature),
                RenderDna::feature_range_text(feature),
            ];
            for key in [
                "gentle_track_source",
                "gentle_track_name",
                "gentle_track_file",
                "gentle_array_dataset",
                "gentle_array_platform",
                "gentle_array_coordinate_system",
                "gentle_array_supported_genome_ids",
                "gentle_array_anchor_genome_id",
                "gentle_array_assembly_check",
                "gentle_array_projection_method",
                "gentle_array_projection_status",
                "gentle_array_native_chromosome",
                "gentle_array_contrast",
                "gentle_array_level",
                "gentle_array_feature_id",
                "gentle_array_value_summary",
                "feature_id",
                "transcript_cluster_id",
                "exon_id",
                "logFC",
                "adj_P_Val",
                "name",
                "label",
                "regulatory_class",
                "note",
                "function",
                "product",
                "gene",
                "id",
                "transcript_id",
                "transcript_variant",
                "variant",
                "repName",
                "repClass",
                "repFamily",
                "rmsk_name",
                "rmsk_class",
                "rmsk_family",
                "repeat_name",
                "repeat_class",
                "repeat_family",
                "rpt_type",
                "rpt_family",
                "mobile_element_type",
            ] {
                for value in feature.qualifier_values(key) {
                    let trimmed = value.trim();
                    if !trimmed.is_empty() {
                        values.push(trimmed.to_string());
                    }
                }
            }
            if RenderDna::is_repeat_feature(feature) {
                values.push("repeat".to_string());
                values.push("rmsk".to_string());
            }
            if RenderDna::is_track_feature(feature) {
                values.push("track".to_string());
                values.push("tracks".to_string());
            }
            if RenderDna::is_array_track_feature(feature) {
                values.push("array".to_string());
                values.push("track".to_string());
            }
            if RenderDna::is_regulatory_feature(feature)
                && !RenderDna::is_track_feature(feature)
                && !RenderDna::is_tfbs_feature(feature)
            {
                values.push("regulatory".to_string());
            }
            values
                .into_iter()
                .map(|value| value.to_ascii_lowercase())
                .collect::<Vec<_>>()
        };
        let source_terms = feature
            .qualifier_values("gentle_track_source")
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();
        let track_terms = feature
            .qualifier_values("gentle_track_name")
            .chain(feature.qualifier_values("gentle_array_platform"))
            .chain(feature.qualifier_values("name"))
            .chain(feature.qualifier_values("label"))
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();
        let file_terms = feature
            .qualifier_values("gentle_track_file")
            .flat_map(|value| {
                let mut values = Vec::new();
                let trimmed = value.trim();
                if !trimmed.is_empty() {
                    values.push(trimmed.to_ascii_lowercase());
                    if let Some(file_name) = Path::new(trimmed).file_name().and_then(|v| v.to_str())
                    {
                        values.push(file_name.to_ascii_lowercase());
                    }
                }
                values
            })
            .collect::<Vec<_>>();
        let note_terms = feature
            .qualifier_values("note")
            .chain(feature.qualifier_values("function"))
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();
        let mut regulatory_terms = feature
            .qualifier_values("regulatory_class")
            .chain(feature.qualifier_values("note"))
            .chain(feature.qualifier_values("function"))
            .chain(feature.qualifier_values("label"))
            .chain(feature.qualifier_values("name"))
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();
        if RenderDna::is_regulatory_feature(feature)
            && !RenderDna::is_track_feature(feature)
            && !RenderDna::is_tfbs_feature(feature)
        {
            regulatory_terms.push("regulatory".to_string());
        }
        let array_terms = feature
            .qualifier_values("gentle_array_dataset")
            .chain(feature.qualifier_values("gentle_array_platform"))
            .chain(feature.qualifier_values("gentle_array_coordinate_system"))
            .chain(feature.qualifier_values("gentle_array_supported_genome_ids"))
            .chain(feature.qualifier_values("gentle_array_anchor_genome_id"))
            .chain(feature.qualifier_values("gentle_array_assembly_check"))
            .chain(feature.qualifier_values("gentle_array_projection_method"))
            .chain(feature.qualifier_values("gentle_array_projection_status"))
            .chain(feature.qualifier_values("gentle_array_native_chromosome"))
            .chain(feature.qualifier_values("gentle_array_contrast"))
            .chain(feature.qualifier_values("gentle_array_level"))
            .chain(feature.qualifier_values("gentle_array_feature_id"))
            .chain(feature.qualifier_values("gentle_array_value_summary"))
            .chain(feature.qualifier_values("feature_id"))
            .chain(feature.qualifier_values("transcript_cluster_id"))
            .chain(feature.qualifier_values("exon_id"))
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();
        let contrast_terms = feature
            .qualifier_values("gentle_array_contrast")
            .chain(feature.qualifier_values("gentle_array_value_summary"))
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();
        let kind_terms = vec![
            feature.kind.to_string().to_ascii_lowercase(),
            kind_label.trim().to_ascii_lowercase(),
        ];
        let mut label_terms = feature
            .qualifier_values("label")
            .chain(feature.qualifier_values("name"))
            .chain(feature.qualifier_values("gene"))
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();
        if !display_label.trim().is_empty() {
            label_terms.push(display_label.trim().to_ascii_lowercase());
        }
        let mut range_terms = vec![RenderDna::feature_range_text(feature).to_ascii_lowercase()];
        if !range_label.trim().is_empty() {
            range_terms.push(range_label.trim().to_ascii_lowercase());
        }
        let mut repeat_terms = Vec::new();
        if let Some(label) = RenderDna::repeat_group_label(feature) {
            repeat_terms.push(label.to_ascii_lowercase());
        }
        for key in [
            "repName",
            "repClass",
            "repFamily",
            "rmsk_name",
            "rmsk_class",
            "rmsk_family",
            "repeat_name",
            "repeat_class",
            "repeat_family",
            "rpt_type",
            "rpt_family",
            "mobile_element_type",
        ] {
            for value in feature.qualifier_values(key) {
                let trimmed = value.trim();
                if !trimmed.is_empty() {
                    repeat_terms.push(trimmed.to_ascii_lowercase());
                }
            }
        }

        terms.iter().all(|(scope, needle)| {
            let search_space: &Vec<String> = match scope.as_deref() {
                Some("source") | Some("type") => &source_terms,
                Some("track") | Some("name") => &track_terms,
                Some("path") | Some("file") => &file_terms,
                Some("note") => &note_terms,
                Some("regulatory") | Some("reg") => &regulatory_terms,
                Some("array") => &array_terms,
                Some("contrast") => &contrast_terms,
                Some("kind") => &kind_terms,
                Some("label") => &label_terms,
                Some("range") => &range_terms,
                Some("repeat") | Some("repeat_class") | Some("rmsk") => &repeat_terms,
                _ => &all_terms,
            };
            search_space.iter().any(|value| value.contains(needle))
        })
    }

    pub(super) fn parse_feature_tree_filter_terms(
        filter_text: &str,
    ) -> Vec<(Option<String>, String)> {
        filter_text
            .split_whitespace()
            .filter_map(|raw_term| {
                let term = raw_term.trim();
                if term.is_empty() {
                    return None;
                }
                if let Some((raw_scope, raw_needle)) = term.split_once(':') {
                    let scope = raw_scope.trim().to_ascii_lowercase();
                    let needle = raw_needle.trim().to_ascii_lowercase();
                    if !scope.is_empty() && !needle.is_empty() {
                        return Some((Some(scope), needle));
                    }
                }
                Some((None, term.to_ascii_lowercase()))
            })
            .collect()
    }

    pub(super) fn feature_tree_filter_help_text() -> &'static str {
        "Free text matches kind/label/range and qualifiers. Scoped terms: kind:mrna label:tp73 range:6128..16430 regulatory:enhancer track:Clariom contrast:AdTAp73alpha-AdGFP source:array path:peaks.bed note:enhancer repeat_class:line"
    }
    pub(super) fn append_filter_term(filter_text: &mut String, term: &str) {
        let normalized = term.trim();
        if normalized.is_empty() {
            return;
        }
        if filter_text
            .split_whitespace()
            .any(|existing| existing.eq_ignore_ascii_case(normalized))
        {
            return;
        }
        if !filter_text.trim().is_empty() {
            filter_text.push(' ');
        }
        filter_text.push_str(normalized);
    }

    pub(super) fn filter_term_present(filter_text: &str, term: &str) -> bool {
        let normalized = term.trim();
        if normalized.is_empty() {
            return false;
        }
        filter_text
            .split_whitespace()
            .any(|existing| existing.eq_ignore_ascii_case(normalized))
    }

    pub(super) fn remove_filter_term(filter_text: &mut String, term: &str) {
        let normalized = term.trim();
        if normalized.is_empty() {
            return;
        }
        let kept = filter_text
            .split_whitespace()
            .filter(|existing| !existing.eq_ignore_ascii_case(normalized))
            .collect::<Vec<_>>();
        *filter_text = kept.join(" ");
    }

    pub(super) fn set_filter_term_enabled(filter_text: &mut String, term: &str, enabled: bool) {
        if enabled {
            Self::append_filter_term(filter_text, term);
        } else {
            Self::remove_filter_term(filter_text, term);
        }
    }

    pub(super) fn apply_anchored_promoter_preset_defaults(&mut self) {
        self.anchored_mode_feature = true;
        self.anchored_feature_kind = "CDS".to_string();
        self.anchored_feature_label.clear();
        self.anchored_feature_boundary_start = true;
        self.anchored_feature_occurrence = "0".to_string();
        self.anchored_direction_upstream = true;
        self.anchored_target_len = "500".to_string();
        self.anchored_tolerance = "100".to_string();
        self.anchored_required_re_sites.clear();
        self.anchored_required_tf_motifs.clear();
        self.anchored_forward_primer.clear();
        self.anchored_reverse_primer.clear();
        self.anchored_output_prefix = "promoter".to_string();
        self.anchored_unique = false;
        self.anchored_max_candidates = "20".to_string();
    }

    pub(super) fn normalize_anchor_output_token(raw: &str) -> String {
        let mut out = String::new();
        let mut previous_underscore = false;
        for ch in raw.chars() {
            if ch.is_ascii_alphanumeric() {
                out.push(ch.to_ascii_lowercase());
                previous_underscore = false;
            } else if !previous_underscore {
                out.push('_');
                previous_underscore = true;
            }
        }
        out.trim_matches('_').to_string()
    }

    pub(super) fn seed_anchored_promoter_from_feature_id(&mut self, feature_id: usize) {
        let feature = {
            self.dna
                .read()
                .ok()
                .and_then(|dna| dna.features().get(feature_id).cloned())
        };
        let Some(feature) = feature else {
            self.op_status =
                format!("Could not seed promoter anchor: feature #{feature_id} not found");
            return;
        };

        self.apply_anchored_promoter_preset_defaults();
        self.anchored_feature_kind = feature.kind.to_string();
        self.anchored_feature_label = Self::feature_tree_first_nonempty_qualifier(
            &feature,
            &["gene", "label", "gene_id", "transcript_id", "name"],
        )
        .or_else(|| {
            let name = RenderDna::feature_name(&feature);
            let trimmed = name.trim();
            (!trimmed.is_empty()).then_some(trimmed.to_string())
        })
        .unwrap_or_default();

        let is_reverse = feature_is_reverse(&feature);
        self.anchored_feature_boundary_start = !is_reverse;
        self.anchored_direction_upstream = !is_reverse;
        let label_token = Self::normalize_anchor_output_token(&self.anchored_feature_label);
        if !label_token.is_empty() {
            self.anchored_output_prefix = format!("promoter_{label_token}");
        }
        self.op_status = format!(
            "Seeded anchored promoter extraction from feature #{} (kind='{}', strand='{}', boundary='{}', direction='{}')",
            feature_id + 1,
            self.anchored_feature_kind,
            if is_reverse { "-" } else { "+" },
            if self.anchored_feature_boundary_start {
                "start"
            } else {
                "end"
            },
            if self.anchored_direction_upstream {
                "upstream"
            } else {
                "downstream"
            }
        );
        self.save_engine_ops_state();
    }

    pub(super) fn current_feature_tree_cache_key(
        &self,
        viewport: Option<(usize, usize)>,
    ) -> FeatureTreeCacheKey {
        let (seq_len, feature_count) = self
            .dna
            .read()
            .map(|dna| (dna.len(), dna.features().len()))
            .unwrap_or((0, 0));
        let (
            show_cds_features,
            show_gene_features,
            show_mrna_features,
            show_repeat_features,
            show_array_features,
            show_contextual_transcript_features,
            show_tfbs,
            tfbs_display_criteria,
            vcf_display_criteria,
            hidden_feature_kinds,
        ) = {
            let display = self.dna_display.read().expect("DNA display lock poisoned");
            (
                display.show_cds_features_effective(),
                display.show_gene_features(),
                display.show_mrna_features(),
                display.show_repeat_features(),
                display.show_array_features(),
                display.show_contextual_transcript_features(),
                display.show_tfbs(),
                display.tfbs_display_criteria(),
                display.vcf_display_criteria(),
                display.hidden_feature_kinds().clone(),
            )
        };
        FeatureTreeCacheKey {
            seq_len,
            feature_count,
            viewport,
            grouping_mode: self.feature_tree_grouping_mode,
            feature_filter_text: self.feature_tree_filter.trim().to_string(),
            show_cds_features,
            show_gene_features,
            show_mrna_features,
            show_repeat_features,
            show_array_features,
            show_contextual_transcript_features,
            show_tfbs,
            tfbs_display_criteria,
            vcf_display_criteria,
            hidden_feature_kinds,
        }
    }

    pub(super) fn build_feature_tree_model(
        &self,
        key: &FeatureTreeCacheKey,
    ) -> FeatureTreeComputedModel {
        let filter_active = !key.feature_filter_text.is_empty();
        let mut filter_total_count = 0usize;
        let mut filter_matched_count = 0usize;
        let typed_features = {
            let dna = self.dna.read().expect("DNA lock poisoned");
            let sequence_length = dna.len();
            dna.features()
                .iter()
                .enumerate()
                .filter_map(|(id, feature)| {
                    if RenderDna::is_source_feature(feature) {
                        return None;
                    }
                    if !RenderDna::feature_passes_kind_filter(
                        feature,
                        key.show_cds_features,
                        key.show_gene_features,
                        key.show_mrna_features,
                    ) {
                        return None;
                    }
                    if !key.show_contextual_transcript_features
                        && RenderDna::is_contextual_transcript_feature(feature)
                    {
                        return None;
                    }
                    if !key.show_repeat_features && RenderDna::is_repeat_feature(feature) {
                        return None;
                    }
                    if !key.show_array_features && RenderDna::is_array_track_feature(feature) {
                        return None;
                    }
                    if key
                        .hidden_feature_kinds
                        .contains(&feature.kind.to_string().to_ascii_uppercase())
                    {
                        return None;
                    }
                    if RenderDna::is_tfbs_feature(feature) {
                        if !key.show_tfbs {
                            return None;
                        }
                        if !RenderDna::tfbs_feature_passes_display_filter(
                            feature,
                            key.tfbs_display_criteria,
                        ) {
                            return None;
                        }
                    }
                    if RenderDna::is_vcf_track_feature(feature)
                        && !RenderDna::vcf_feature_passes_display_filter(
                            feature,
                            &key.vcf_display_criteria,
                        )
                    {
                        return None;
                    }
                    let (from, to) = feature.location.find_bounds().ok()?;
                    if from < 0 || to < 0 {
                        return None;
                    }
                    let visible_in_view = match key.viewport {
                        Some((start, end)) => Self::feature_overlaps_linear_viewport(
                            feature,
                            sequence_length,
                            start,
                            end,
                        ),
                        None => true,
                    };
                    let kind_upper = feature.kind.to_string().to_ascii_uppercase();
                    let is_repeat = RenderDna::is_repeat_feature(feature);
                    let is_track = RenderDna::is_track_feature(feature);
                    let is_array_track = RenderDna::is_array_track_feature(feature);
                    let is_tfbs = RenderDna::is_tfbs_feature(feature);
                    let kind_label = if RenderDna::is_repeat_feature(feature) {
                        feature.kind.to_string()
                    } else if is_array_track {
                        "array".to_string()
                    } else if is_track {
                        "tracks".to_string()
                    } else {
                        feature.kind.to_string()
                    };
                    let base_label = {
                        let name = RenderDna::feature_name(feature);
                        if name.trim().is_empty() {
                            format!("{} #{}", feature.kind, id + 1)
                        } else {
                            name
                        }
                    };
                    let range_label = RenderDna::feature_range_text(feature);
                    let full_display_label =
                        Self::feature_tree_display_label(feature, &base_label, &range_label);
                    let regulatory_grouping =
                        Self::derive_regulatory_feature_grouping(feature, &full_display_label);
                    let display_label = regulatory_grouping
                        .as_ref()
                        .map(|grouping| grouping.display_label.clone())
                        .unwrap_or_else(|| full_display_label.clone());
                    if filter_active {
                        filter_total_count += 1;
                        if !Self::feature_tree_matches_filter(
                            feature,
                            &key.feature_filter_text,
                            &kind_label,
                            &full_display_label,
                            &range_label,
                        ) {
                            return None;
                        }
                        filter_matched_count += 1;
                    }
                    let subgroup_label =
                        Self::feature_tree_subgroup_label(feature, &base_label, key.grouping_mode);
                    let subgroup_key = subgroup_label
                        .as_ref()
                        .map(|label| label.trim().to_ascii_uppercase());
                    let can_seed_promoter_anchor = kind_label.eq_ignore_ascii_case("mrna")
                        || kind_label.eq_ignore_ascii_case("transcript");
                    let disable_grouping = feature.kind.to_string().eq_ignore_ascii_case("GENE");
                    Some((
                        kind_label,
                        FeatureTreeEntry {
                            id,
                            feature_label_full: full_display_label,
                            feature_label: display_label,
                            range_label,
                            subgroup_key,
                            subgroup_label,
                            prefer_grouped_label: kind_upper.contains("RNA"),
                            show_range_inline_when_ungrouped: matches!(
                                kind_upper.as_str(),
                                "MRNA" | "GENE"
                            ),
                            visible_in_view,
                            is_regulatory: RenderDna::is_regulatory_feature(feature),
                            is_repeat,
                            is_track,
                            is_array_track,
                            is_tfbs,
                            disable_grouping,
                            supports_splicing_expert: Self::feature_kind_supports_splicing_expert(
                                kind_upper.as_str(),
                            ),
                            supports_variant_followup: Self::feature_kind_supports_variant_followup(
                                kind_upper.as_str(),
                            ),
                            can_seed_promoter_anchor,
                            regulatory_primary_group_key: regulatory_grouping
                                .as_ref()
                                .map(|grouping| grouping.primary_key.clone()),
                            regulatory_primary_group_label: regulatory_grouping
                                .as_ref()
                                .map(|grouping| grouping.primary_label.clone()),
                            regulatory_secondary_group_key: regulatory_grouping
                                .as_ref()
                                .and_then(|grouping| grouping.secondary_key.clone()),
                            regulatory_secondary_group_label: regulatory_grouping
                                .as_ref()
                                .and_then(|grouping| grouping.secondary_label.clone()),
                        },
                    ))
                })
                .collect::<Vec<_>>()
        };

        let mut grouped_features: HashMap<String, Vec<FeatureTreeEntry>> = HashMap::new();
        for (kind, entry) in typed_features {
            grouped_features.entry(kind).or_default().push(entry);
        }

        let mut group_keys = grouped_features.keys().cloned().collect::<Vec<_>>();
        group_keys.sort();
        let mut groups = Vec::with_capacity(group_keys.len());
        for kind in group_keys {
            let Some(entries) = grouped_features.remove(&kind) else {
                continue;
            };
            let visible_count = entries.iter().filter(|entry| entry.visible_in_view).count();
            let mut subgroup_cardinality: HashMap<String, usize> = HashMap::new();
            for entry in &entries {
                if let Some(subgroup_key) = &entry.subgroup_key {
                    *subgroup_cardinality
                        .entry(subgroup_key.clone())
                        .or_insert(0) += 1;
                }
            }

            let mut grouped_entry_map: HashMap<String, (String, Vec<usize>)> = HashMap::new();
            let mut regulatory_primary_map: HashMap<String, (String, Vec<usize>)> = HashMap::new();
            let mut ungrouped_entry_indices: Vec<usize> = Vec::new();
            for (index, entry) in entries.iter().enumerate() {
                if matches!(key.grouping_mode, FeatureTreeGroupingMode::Off)
                    || entry.disable_grouping
                {
                    ungrouped_entry_indices.push(index);
                    continue;
                }
                if let (Some(primary_key), Some(primary_label)) = (
                    &entry.regulatory_primary_group_key,
                    &entry.regulatory_primary_group_label,
                ) {
                    regulatory_primary_map
                        .entry(primary_key.clone())
                        .or_insert_with(|| (primary_label.clone(), Vec::new()))
                        .1
                        .push(index);
                    continue;
                }
                let subgroup_count = entry
                    .subgroup_key
                    .as_ref()
                    .and_then(|subgroup_key| subgroup_cardinality.get(subgroup_key))
                    .copied()
                    .unwrap_or(0);
                let should_group =
                    Self::feature_tree_should_group(key.grouping_mode, subgroup_count);
                if let (Some(subgroup_key), Some(subgroup_label)) =
                    (&entry.subgroup_key, &entry.subgroup_label)
                {
                    if should_group {
                        grouped_entry_map
                            .entry(subgroup_key.clone())
                            .or_insert_with(|| (subgroup_label.clone(), Vec::new()))
                            .1
                            .push(index);
                    } else {
                        ungrouped_entry_indices.push(index);
                    }
                } else {
                    ungrouped_entry_indices.push(index);
                }
            }

            let mut grouped_entry_keys = grouped_entry_map.keys().cloned().collect::<Vec<_>>();
            grouped_entry_keys.sort_by(|left, right| {
                let left_label = grouped_entry_map
                    .get(left)
                    .map(|(label, _)| label.as_str())
                    .unwrap_or("");
                let right_label = grouped_entry_map
                    .get(right)
                    .map(|(label, _)| label.as_str())
                    .unwrap_or("");
                left_label.cmp(right_label).then_with(|| left.cmp(right))
            });
            let grouped_entries = grouped_entry_keys
                .into_iter()
                .filter_map(|group_key| {
                    let (label, entry_indices) = grouped_entry_map.remove(&group_key)?;
                    let visible_count = entry_indices
                        .iter()
                        .filter(|index| entries[**index].visible_in_view)
                        .count();
                    Some(FeatureTreeNamedGroup {
                        key: group_key,
                        label,
                        entry_indices,
                        visible_count,
                    })
                })
                .collect::<Vec<_>>();

            let mut regulatory_primary_keys =
                regulatory_primary_map.keys().cloned().collect::<Vec<_>>();
            regulatory_primary_keys.sort_by(|left, right| {
                let left_label = regulatory_primary_map
                    .get(left)
                    .map(|(label, _)| label.as_str())
                    .unwrap_or("");
                let right_label = regulatory_primary_map
                    .get(right)
                    .map(|(label, _)| label.as_str())
                    .unwrap_or("");
                left_label
                    .to_ascii_lowercase()
                    .cmp(&right_label.to_ascii_lowercase())
                    .then_with(|| left.cmp(right))
            });
            let regulatory_primary_groups = regulatory_primary_keys
                .into_iter()
                .filter_map(|primary_key| {
                    let (label, entry_indices) = regulatory_primary_map.remove(&primary_key)?;
                    let visible_count = entry_indices
                        .iter()
                        .filter(|index| entries[**index].visible_in_view)
                        .count();
                    let mut secondary_map: HashMap<String, (String, Vec<usize>)> = HashMap::new();
                    let mut secondary_ungrouped = Vec::new();
                    for index in entry_indices.iter().copied() {
                        let entry = &entries[index];
                        if let (Some(secondary_key), Some(secondary_label)) = (
                            &entry.regulatory_secondary_group_key,
                            &entry.regulatory_secondary_group_label,
                        ) {
                            secondary_map
                                .entry(secondary_key.clone())
                                .or_insert_with(|| (secondary_label.clone(), Vec::new()))
                                .1
                                .push(index);
                        } else {
                            secondary_ungrouped.push(index);
                        }
                    }
                    let mut secondary_keys = secondary_map.keys().cloned().collect::<Vec<_>>();
                    secondary_keys.sort_by(|left, right| {
                        let left_label = secondary_map
                            .get(left)
                            .map(|(label, _)| label.as_str())
                            .unwrap_or("");
                        let right_label = secondary_map
                            .get(right)
                            .map(|(label, _)| label.as_str())
                            .unwrap_or("");
                        left_label
                            .to_ascii_lowercase()
                            .cmp(&right_label.to_ascii_lowercase())
                            .then_with(|| left.cmp(right))
                    });
                    let secondary_groups = secondary_keys
                        .into_iter()
                        .filter_map(|secondary_key| {
                            let (secondary_label, secondary_entry_indices) =
                                secondary_map.remove(&secondary_key)?;
                            let visible_count = secondary_entry_indices
                                .iter()
                                .filter(|index| entries[**index].visible_in_view)
                                .count();
                            Some(FeatureTreeNamedGroup {
                                key: secondary_key,
                                label: secondary_label,
                                entry_indices: secondary_entry_indices,
                                visible_count,
                            })
                        })
                        .collect::<Vec<_>>();
                    Some(FeatureTreeRegulatoryPrimaryGroup {
                        key: primary_key,
                        label,
                        entry_indices,
                        visible_count,
                        secondary_groups,
                        ungrouped_entry_indices: secondary_ungrouped,
                    })
                })
                .collect::<Vec<_>>();

            groups.push(FeatureTreeKindGroup {
                kind,
                entries,
                visible_count,
                grouped_entries,
                regulatory_primary_groups,
                ungrouped_entry_indices,
            });
        }

        FeatureTreeComputedModel {
            filter_total_count,
            filter_matched_count,
            groups,
        }
    }

    pub(super) fn ensure_feature_tree_cache_current(&mut self, viewport: Option<(usize, usize)>) {
        let next_key = self.current_feature_tree_cache_key(viewport);
        let is_current = self
            .feature_tree_cache
            .as_ref()
            .map(|cache| cache.key == next_key)
            .unwrap_or(false);
        if is_current {
            return;
        }
        let model = self.build_feature_tree_model(&next_key);
        self.feature_tree_cache = Some(FeatureTreeCache {
            key: next_key,
            model,
        });
    }

    pub(super) fn render_features(&mut self, ui: &mut egui::Ui) {
        ui.heading(
            self.dna
                .read()
                .expect("DNA lock poisoned")
                .name()
                .as_ref()
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
                .unwrap_or_else(|| {
                    self.seq_id
                        .clone()
                        .unwrap_or_else(|| "<Unnamed DNA sequence>".to_string())
                }),
        );
        let multi_select_count = self.multi_selected_feature_ids.len();
        if multi_select_count > 1 {
            ui.horizontal_wrapped(|ui| {
                ui.label(
                    egui::RichText::new(format!("Multi-select active ({multi_select_count})"))
                        .size(10.5)
                        .strong()
                        .color(egui::Color32::from_rgb(7, 56, 96))
                        .background_color(egui::Color32::from_rgb(218, 234, 248)),
                )
                .on_hover_text("Multiple features are selected and forced to use external labels");
                if ui
                    .small_button("Clear multi-select")
                    .on_hover_text(
                        "Keep the currently focused feature selected and remove other selections",
                    )
                    .clicked()
                {
                    self.clear_multi_select_keep_primary();
                }
            });
        }
        let seq_key = self.panel_scope_key();
        let viewport = self.active_linear_viewport_range();
        let viewport_limited = viewport.is_some();
        ui.horizontal(|ui| {
            ui.label(Self::tr("sequence.grouping"));
            let before = self.feature_tree_grouping_mode;
            egui::ComboBox::from_id_salt(format!("feature_tree_grouping_mode_{seq_key}"))
                .selected_text(self.feature_tree_grouping_mode.label())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.feature_tree_grouping_mode,
                        FeatureTreeGroupingMode::Off,
                        FeatureTreeGroupingMode::Off.label(),
                    );
                    ui.selectable_value(
                        &mut self.feature_tree_grouping_mode,
                        FeatureTreeGroupingMode::Auto,
                        FeatureTreeGroupingMode::Auto.label(),
                    );
                    ui.selectable_value(
                        &mut self.feature_tree_grouping_mode,
                        FeatureTreeGroupingMode::Always,
                        FeatureTreeGroupingMode::Always.label(),
                    );
                });
            if before != self.feature_tree_grouping_mode {
                self.save_engine_ops_state();
            }
        })
        .response
        .on_hover_text(
            "Off: flat list per kind. Auto: subgroup repeated labels only. Always: subgroup all labels when available.",
        );
        ui.horizontal(|ui| {
            ui.label(Self::tr("sequence.filter"))
                .on_hover_text(Self::feature_tree_filter_help_text());
            let response = ui
                .add(
                    egui::TextEdit::singleline(&mut self.feature_tree_filter)
                        .hint_text("kind:mrna label:tp73 range:6128..16430")
                        .desired_width(220.0),
                )
                .on_hover_text(Self::feature_tree_filter_help_text());
            if response.changed() {
                self.pending_feature_tree_scroll_to = None;
                self.save_engine_ops_state();
            }
            if ui
                .add_enabled(
                    !self.feature_tree_filter.trim().is_empty(),
                    egui::Button::new(Self::tr("button.clear")),
                )
                .clicked()
            {
                self.feature_tree_filter.clear();
                self.pending_feature_tree_scroll_to = None;
                self.save_engine_ops_state();
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.small("Focus:");
            for preset in [
                FeatureTreeFocusPreset::Cloning,
                FeatureTreeFocusPreset::Regulatory,
                FeatureTreeFocusPreset::Repeats,
                FeatureTreeFocusPreset::Tracks,
            ] {
                let active = self.feature_tree_filter.trim() == preset.filter_text()
                    && self.feature_tree_grouping_mode == preset.grouping_mode();
                if ui
                    .selectable_label(active, preset.label())
                    .on_hover_text(preset.hover_text())
                    .clicked()
                {
                    self.apply_feature_tree_focus_preset(preset);
                }
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.small(Self::tr("sequence.presets"));
            let mut presets_changed = false;
            for preset in [
                "kind:tracks",
                "kind:array",
                "source:array",
                "source:bed",
                "source:vcf",
                "track:Clariom",
                "path:.bed",
                "note:enhancer",
                "label:tp73",
            ] {
                let active = Self::filter_term_present(&self.feature_tree_filter, preset);
                let response = ui
                    .selectable_label(active, preset)
                    .on_hover_text(if active {
                        "Disable this preset filter term"
                    } else {
                        "Enable this preset filter term"
                    });
                if response.clicked() {
                    Self::set_filter_term_enabled(&mut self.feature_tree_filter, preset, !active);
                    presets_changed = true;
                }
            }
            if presets_changed {
                self.pending_feature_tree_scroll_to = None;
                self.save_engine_ops_state();
            }
        });
        self.ensure_feature_tree_cache_current(viewport);
        if let Some(model) = self.feature_tree_cache.as_ref().map(|cache| &cache.model) {
            let layer_labels = Self::feature_tree_layer_summary_labels(model, viewport_limited);
            if !layer_labels.is_empty() {
                ui.horizontal_wrapped(|ui| {
                    ui.small("Layers:");
                    for label in layer_labels {
                        ui.label(egui::RichText::new(label).size(10.5))
                            .on_hover_text(
                                "Counts reflect the current feature-tree filter and linear view span",
                            );
                    }
                });
            }
        }
        if viewport_limited {
            ui.label(
                egui::RichText::new("Group counts show visible/total for current linear view span")
                    .size(10.5)
                    .italics(),
            );
        }

        let selected_feature_ids = self.multi_selected_feature_ids.clone();
        let feature_details_font_size = self.feature_details_font_size();
        let can_fit_linear = !self.is_circular();
        let mut clicked_feature: Option<(usize, bool)> = None;
        let mut fit_feature: Option<usize> = None;
        let mut promoter_anchor_seed_feature: Option<usize> = None;
        let mut open_variant_followup_feature: Option<usize> = None;
        let mut open_splicing_feature: Option<usize> = None;
        let mut open_rna_read_mapping_feature: Option<usize> = None;
        let mut open_dotplot_feature: Option<usize> = None;
        let mut copy_feature_payload: Option<(usize, FeatureCopyPayloadKind)> = None;
        let mut focus_matching_array_feature: Option<usize> = None;
        let feature_font_size = feature_details_font_size;
        let kind_font_size = feature_font_size + 1.0;
        {
            let Some(model) = self.feature_tree_cache.as_ref().map(|cache| &cache.model) else {
                return;
            };
            if !self.feature_tree_filter.trim().is_empty() {
                ui.small(format!(
                    "Filter matches: {} / {} features",
                    model.filter_matched_count, model.filter_total_count
                ));
            }
            for group in &model.groups {
                let has_selected = group
                    .entries
                    .iter()
                    .any(|entry| selected_feature_ids.contains(&entry.id));
                egui::CollapsingHeader::new(
                    egui::RichText::new(Self::format_feature_tree_count_label(
                        group.kind.as_str(),
                        group.visible_count,
                        group.entries.len(),
                        viewport_limited,
                    ))
                    .size(kind_font_size)
                    .strong(),
                )
                .id_salt(format!("feature_kind_{seq_key}_{}", group.kind))
                .open(if has_selected { Some(true) } else { None })
                .show(ui, |ui| {
                    let mut render_entry =
                        |ui: &mut egui::Ui,
                         entry: &FeatureTreeEntry,
                         grouped_entry: bool,
                         show_range_only_when_grouped: bool| {
                            let selected = selected_feature_ids.contains(&entry.id);
                            let button_label = if grouped_entry
                                && show_range_only_when_grouped
                                && !entry.prefer_grouped_label
                                && !entry.range_label.is_empty()
                            {
                                entry.range_label.clone()
                            } else if entry.show_range_inline_when_ungrouped
                                && !grouped_entry
                                && !entry.range_label.is_empty()
                            {
                                format!("{} [{}]", entry.feature_label, entry.range_label)
                            } else {
                                entry.feature_label.clone()
                            };
                            let mut button_text =
                                egui::RichText::new(button_label).size(feature_font_size);
                            if viewport_limited && !entry.visible_in_view {
                                button_text = button_text.color(egui::Color32::from_gray(120));
                            }
                            ui.horizontal(|ui| {
                                let button = egui::Button::new(button_text).selected(selected);
                                let mut response = ui.add(button);
                                let tail_width = ui.available_width().max(0.0);
                                if tail_width > 1.0 {
                                    let tail_height =
                                        response.rect.height().max(ui.spacing().interact_size.y);
                                    response = response.union(ui.allocate_response(
                                        egui::vec2(tail_width, tail_height),
                                        egui::Sense::click(),
                                    ));
                                }
                                let mut hover_lines: Vec<String> = Vec::new();
                                let feature_label_with_range = if entry.range_label.is_empty() {
                                    entry.feature_label_full.clone()
                                } else {
                                    format!("{} ({})", entry.feature_label_full, entry.range_label)
                                };
                                if grouped_entry && show_range_only_when_grouped {
                                    hover_lines.push(feature_label_with_range);
                                } else if entry.feature_label_full != entry.feature_label {
                                    hover_lines.push(feature_label_with_range);
                                }
                                if entry.is_regulatory
                                    && entry.feature_label != entry.feature_label_full
                                {
                                    let hover_text = if entry.range_label.is_empty() {
                                        entry.feature_label.clone()
                                    } else {
                                        format!("{} ({})", entry.feature_label, entry.range_label)
                                    };
                                    hover_lines.push(format!("Display: {hover_text}"));
                                }
                                if viewport_limited && !entry.visible_in_view {
                                    hover_lines.push("Outside current linear view span".to_string());
                                }
                                if !hover_lines.is_empty() {
                                    response = response.on_hover_text(hover_lines.join("\n"));
                                }
                                if selected && self.pending_feature_tree_scroll_to == Some(entry.id) {
                                    response.scroll_to_me(Some(egui::Align::Center));
                                    self.pending_feature_tree_scroll_to = None;
                                }
                                if response.clicked() {
                                    let additive = ui.input(|i| i.modifiers.command);
                                    clicked_feature = Some((entry.id, additive));
                                }
                                if response.double_clicked() {
                                    clicked_feature = Some((entry.id, false));
                                    open_splicing_feature = Some(entry.id);
                                }
                                response.context_menu(|ui| {
                                    if ui
                                        .button("Focus feature (current zoom)")
                                        .on_hover_text("Center selected feature without changing zoom")
                                        .clicked()
                                    {
                                        clicked_feature = Some((entry.id, false));
                                        ui.close();
                                    }
                                    let fit_response = ui.add_enabled(
                                        can_fit_linear,
                                        egui::Button::new("Fit feature in view"),
                                    );
                                    let fit_response = if can_fit_linear {
                                        fit_response.on_hover_text(
                                            "Adjust linear viewport so this feature is fully visible",
                                        )
                                    } else {
                                        fit_response.on_hover_text(
                                            "Available in linear map mode only",
                                        )
                                    };
                                    if fit_response.clicked() {
                                        clicked_feature = Some((entry.id, false));
                                        fit_feature = Some(entry.id);
                                        ui.close();
                                    }
                                    ui.separator();
                                    if entry.is_array_track {
                                        if ui
                                            .button("Copy array value table")
                                            .on_hover_text(
                                                "Copy the projected array values and metadata for this probeset",
                                            )
                                            .clicked()
                                        {
                                            copy_feature_payload =
                                                Some((entry.id, FeatureCopyPayloadKind::PopupText));
                                            ui.close();
                                            return;
                                        }
                                        if ui
                                            .button("Focus matching probesets")
                                            .on_hover_text(
                                                "Select matching probeset intervals across all contrast lanes",
                                            )
                                            .clicked()
                                        {
                                            clicked_feature = Some((entry.id, false));
                                            focus_matching_array_feature = Some(entry.id);
                                            ui.close();
                                            return;
                                        }
                                        ui.separator();
                                    }
                                    for kind in [
                                        FeatureCopyPayloadKind::Identifier,
                                        FeatureCopyPayloadKind::Description,
                                        FeatureCopyPayloadKind::PopupText,
                                    ] {
                                        if ui
                                            .button(kind.menu_label())
                                            .on_hover_text(kind.hover_text())
                                            .clicked()
                                        {
                                            copy_feature_payload = Some((entry.id, kind));
                                            ui.close();
                                            return;
                                        }
                                    }
                                    ui.separator();
                                    let promoter_response = ui.add_enabled(
                                        entry.can_seed_promoter_anchor,
                                        egui::Button::new("Use as promoter anchor (Engine Ops)"),
                                    );
                                    let promoter_response = if entry.can_seed_promoter_anchor {
                                        promoter_response.on_hover_text(
                                            "Seed anchored promoter extraction from this mRNA/transcript boundary (strand-aware)",
                                        )
                                    } else {
                                        promoter_response.on_hover_text(
                                            "Available for mRNA/transcript features",
                                        )
                                    };
                                    if promoter_response.clicked() {
                                        clicked_feature = Some((entry.id, false));
                                        promoter_anchor_seed_feature = Some(entry.id);
                                        ui.close();
                                    }
                                    let variant_followup_response = ui.add_enabled(
                                        entry.supports_variant_followup,
                                        egui::Button::new("Open Promoter Design"),
                                    );
                                    let variant_followup_response =
                                        variant_followup_response.on_hover_text(
                                            "Open the dedicated Promoter design window for this promoter-relevant feature",
                                        );
                                    if variant_followup_response.clicked() {
                                        clicked_feature = Some((entry.id, false));
                                        open_variant_followup_feature = Some(entry.id);
                                        ui.close();
                                    }
                                    let splicing_response = ui.add_enabled(
                                        entry.supports_splicing_expert,
                                        egui::Button::new("Open Splicing Window"),
                                    );
                                    let splicing_response = splicing_response.on_hover_text(
                                        "Open the dedicated Splicing Expert window for this feature",
                                    );
                                    if splicing_response.clicked() {
                                        clicked_feature = Some((entry.id, false));
                                        open_splicing_feature = Some(entry.id);
                                        ui.close();
                                    }
                                    let mapping_response = ui.add_enabled(
                                        entry.supports_splicing_expert,
                                        egui::Button::new("Open RNA-read Mapping"),
                                    );
                                    let mapping_response = mapping_response.on_hover_text(
                                        "Open the dedicated RNA-read Mapping workspace for this feature's splicing locus",
                                    );
                                    if mapping_response.clicked() {
                                        clicked_feature = Some((entry.id, false));
                                        open_rna_read_mapping_feature = Some(entry.id);
                                        ui.close();
                                    }
                                    let dotplot_response = ui.add_enabled(
                                        entry.supports_splicing_expert,
                                        egui::Button::new("Derive + Dotplot"),
                                    );
                                    let dotplot_response = dotplot_response.on_hover_text(
                                        "Derive the selected transcript sequence and open the transcript-vs-genomic dotplot workspace",
                                    );
                                    if dotplot_response.clicked() {
                                        clicked_feature = Some((entry.id, false));
                                        open_dotplot_feature = Some(entry.id);
                                        ui.close();
                                    }
                                });
                            });
                        };

                    for primary_group in &group.regulatory_primary_groups {
                        let primary_has_selected = primary_group
                            .entry_indices
                            .iter()
                            .any(|index| selected_feature_ids.contains(&group.entries[*index].id));
                        let primary_heading = Self::format_feature_tree_count_label(
                            &primary_group.label,
                            primary_group.visible_count,
                            primary_group.entry_indices.len(),
                            viewport_limited,
                        );
                        egui::CollapsingHeader::new(
                            egui::RichText::new(primary_heading)
                                .size(feature_font_size)
                                .strong(),
                        )
                        .id_salt(format!(
                            "feature_kind_reg_primary_{seq_key}_{}_{}",
                            group.kind, primary_group.key
                        ))
                        .open(if primary_has_selected { Some(true) } else { None })
                        .show(ui, |ui| {
                            for secondary_group in &primary_group.secondary_groups {
                                let secondary_has_selected =
                                    secondary_group.entry_indices.iter().any(|index| {
                                        selected_feature_ids
                                            .contains(&group.entries[*index].id)
                                    });
                                let secondary_heading = Self::format_feature_tree_count_label(
                                    &secondary_group.label,
                                    secondary_group.visible_count,
                                    secondary_group.entry_indices.len(),
                                    viewport_limited,
                                );
                                egui::CollapsingHeader::new(
                                    egui::RichText::new(secondary_heading)
                                        .size(feature_font_size)
                                        .strong(),
                                )
                                .id_salt(format!(
                                    "feature_kind_reg_secondary_{seq_key}_{}_{}_{}",
                                    group.kind,
                                    primary_group.key,
                                    secondary_group.key
                                ))
                                .open(if secondary_has_selected {
                                    Some(true)
                                } else {
                                    None
                                })
                                .show(ui, |ui| {
                                    for index in &secondary_group.entry_indices {
                                        render_entry(
                                            ui,
                                            &group.entries[*index],
                                            true,
                                            false,
                                        );
                                    }
                                });
                            }

                            for index in &primary_group.ungrouped_entry_indices {
                                render_entry(
                                    ui,
                                    &group.entries[*index],
                                    false,
                                    false,
                                );
                            }
                        });
                    }

                    for subgroup in &group.grouped_entries {
                        let subgroup_has_selected = subgroup
                            .entry_indices
                            .iter()
                            .any(|index| selected_feature_ids.contains(&group.entries[*index].id));
                        let subgroup_heading = Self::format_feature_tree_count_label(
                            &subgroup.label,
                            subgroup.visible_count,
                            subgroup.entry_indices.len(),
                            viewport_limited,
                        );
                        egui::CollapsingHeader::new(
                            egui::RichText::new(subgroup_heading)
                                .size(feature_font_size)
                                .strong(),
                        )
                        .id_salt(format!(
                            "feature_kind_group_{seq_key}_{}_{}",
                            group.kind, subgroup.key
                        ))
                        .open(if subgroup_has_selected { Some(true) } else { None })
                        .show(ui, |ui| {
                            for index in &subgroup.entry_indices {
                                render_entry(
                                    ui,
                                    &group.entries[*index],
                                    true,
                                    true,
                                );
                            }
                        });
                    }

                    for index in &group.ungrouped_entry_indices {
                        render_entry(
                            ui,
                            &group.entries[*index],
                            false,
                            false,
                        );
                    }
                });
            }
        }
        if let Some((feature_id, kind)) = copy_feature_payload {
            self.copy_feature_payload_to_clipboard(ui.ctx(), feature_id, kind);
        }
        if let Some((id, additive)) = clicked_feature {
            if additive {
                self.toggle_feature_multi_select(id);
            } else {
                self.focus_feature(id);
            }
        }
        if let Some(id) = fit_feature {
            self.fit_feature_in_linear_view(id);
        }
        if let Some(feature_id) = promoter_anchor_seed_feature {
            self.seed_anchored_promoter_from_feature_id(feature_id);
        }
        if let Some(feature_id) = open_variant_followup_feature {
            self.open_variant_followup_for_feature(feature_id, "feature tree action");
        }
        if let Some(feature_id) = open_splicing_feature {
            self.open_splicing_expert_for_feature(feature_id, "feature tree action");
        }
        if let Some(feature_id) = open_rna_read_mapping_feature {
            self.open_rna_read_mapping_for_feature(feature_id, "feature tree action");
        }
        if let Some(feature_id) = open_dotplot_feature {
            self.open_dotplot_for_feature(feature_id, "feature tree action");
        }
        if let Some(feature_id) = focus_matching_array_feature {
            self.focus_matching_array_features(feature_id);
        }
    }
}
