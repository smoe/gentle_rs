use crate::{
    dna_display::{DnaDisplay, TfbsDisplayCriteria},
    dna_sequence::DNAsequence,
    render_dna_circular::RenderDnaCircular,
    render_dna_linear::RenderDnaLinear,
    restriction_enzyme::RestrictionEnzymeKey,
};
use eframe::egui::{self, Color32, PointerState, Rect, Response, Sense, Ui, Vec2, Widget};
use gb_io::seq::Feature;
use std::{
    fmt::Debug,
    panic::{catch_unwind, AssertUnwindSafe},
    sync::{Arc, RwLock},
};

#[derive(Debug, Clone, PartialEq)]
pub struct RestrictionEnzymePosition {
    pub area: Rect,
    pub key: RestrictionEnzymeKey,
}

impl RestrictionEnzymePosition {
    pub fn key(&self) -> &RestrictionEnzymeKey {
        &self.key
    }
    pub fn area(&self) -> Rect {
        self.area
    }
}

#[derive(Debug, Clone)]
pub enum RenderDna {
    Circular(Arc<RwLock<RenderDnaCircular>>),
    Linear(Arc<RwLock<RenderDnaLinear>>),
}

impl RenderDna {
    pub fn new(dna: Arc<RwLock<DNAsequence>>, display: Arc<RwLock<DnaDisplay>>) -> Self {
        let is_circular = dna.read().map(|d| d.is_circular()).unwrap_or(true);
        match is_circular {
            true => {
                RenderDna::Circular(Arc::new(RwLock::new(RenderDnaCircular::new(dna, display))))
            }
            false => RenderDna::Linear(Arc::new(RwLock::new(RenderDnaLinear::new(dna, display)))),
        }
    }

    pub fn is_circular(&self) -> bool {
        match self {
            RenderDna::Circular(_) => true,
            RenderDna::Linear(_) => false,
        }
    }

    pub fn on_click(&self, pointer_state: PointerState) {
        match self {
            RenderDna::Circular(renderer) => {
                if let Ok(mut renderer) = renderer.write() {
                    renderer.on_click(pointer_state);
                }
            }
            RenderDna::Linear(renderer) => {
                if let Ok(mut renderer) = renderer.write() {
                    renderer.on_click(pointer_state);
                }
            }
        }
    }

    pub fn on_hover(&self, pointer_state: PointerState) {
        match self {
            RenderDna::Circular(renderer) => {
                if let Ok(mut renderer) = renderer.write() {
                    renderer.on_hover(pointer_state);
                }
            }
            RenderDna::Linear(renderer) => {
                if let Ok(mut renderer) = renderer.write() {
                    renderer.on_hover(pointer_state);
                }
            }
        }
    }

    pub fn on_double_click(&self, pointer_state: PointerState) {
        match self {
            RenderDna::Circular(renderer) => {
                if let Ok(mut renderer) = renderer.write() {
                    renderer.on_double_click(pointer_state);
                }
            }
            RenderDna::Linear(renderer) => {
                if let Ok(mut renderer) = renderer.write() {
                    renderer.on_double_click(pointer_state);
                }
            }
        }
    }

    pub fn get_selected_feature_id(&self) -> Option<usize> {
        match self {
            RenderDna::Circular(renderer) => renderer
                .read()
                .ok()
                .and_then(|r| r.selected_feature_number()),
            RenderDna::Linear(renderer) => renderer
                .read()
                .ok()
                .and_then(|r| r.selected_feature_number()),
        }
    }

    pub fn get_hovered_feature_id(&self) -> Option<usize> {
        match self {
            RenderDna::Circular(renderer) => renderer
                .read()
                .ok()
                .and_then(|r| r.hovered_feature_number()),
            RenderDna::Linear(renderer) => renderer
                .read()
                .ok()
                .and_then(|r| r.hovered_feature_number()),
        }
    }

    pub fn select_feature(&self, feature_number: Option<usize>) {
        match self {
            RenderDna::Circular(renderer) => {
                if let Ok(mut renderer) = renderer.write() {
                    renderer.select_feature(feature_number);
                }
            }
            RenderDna::Linear(renderer) => {
                if let Ok(mut renderer) = renderer.write() {
                    renderer.select_feature(feature_number);
                }
            }
        }
    }

    fn render(&self, ui: &mut egui::Ui, area: Rect) {
        let result = catch_unwind(AssertUnwindSafe(|| match self {
            RenderDna::Circular(renderer) => {
                if let Ok(mut renderer) = renderer.write() {
                    renderer.render(ui, area);
                }
            }
            RenderDna::Linear(renderer) => {
                if let Ok(mut renderer) = renderer.write() {
                    renderer.render(ui, area);
                }
            }
        }));
        if result.is_err() {
            eprintln!("W RenderDna: recovered from panic while rendering map");
        }
    }

    pub fn is_feature_pointy(feature: &Feature) -> bool {
        matches!(
            feature.kind.to_string().to_ascii_uppercase().as_str(),
            "CDS" | "GENE"
        )
    }

    pub fn feature_color(feature: &Feature) -> Color32 {
        match feature.kind.to_string().to_ascii_uppercase().as_str() {
            "CDS" => Color32::RED,
            "GENE" => Color32::BLUE,
            "MRNA" => Color32::from_rgb(180, 100, 10),
            "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND" => Color32::from_rgb(35, 120, 35),
            "TRACK" => Color32::from_rgb(85, 85, 85),
            _ => Color32::GRAY,
        }
    }

    pub fn is_cds_feature(feature: &Feature) -> bool {
        feature.kind.to_string().to_ascii_uppercase() == "CDS"
    }

    pub fn is_gene_feature(feature: &Feature) -> bool {
        feature.kind.to_string().to_ascii_uppercase() == "GENE"
    }

    pub fn is_mrna_feature(feature: &Feature) -> bool {
        feature.kind.to_string().to_ascii_uppercase() == "MRNA"
    }

    pub fn feature_passes_kind_filter(
        feature: &Feature,
        show_cds_features: bool,
        show_gene_features: bool,
        show_mrna_features: bool,
    ) -> bool {
        if Self::is_cds_feature(feature) && !show_cds_features {
            return false;
        }
        if Self::is_gene_feature(feature) && !show_gene_features {
            return false;
        }
        if Self::is_mrna_feature(feature) && !show_mrna_features {
            return false;
        }
        true
    }

    pub fn is_tfbs_feature(feature: &Feature) -> bool {
        matches!(
            feature.kind.to_string().to_ascii_uppercase().as_str(),
            "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND"
        )
    }

    pub fn is_track_feature(feature: &Feature) -> bool {
        feature
            .qualifier_values("gentle_generated".into())
            .any(|value| {
                matches!(
                    value.trim().to_ascii_lowercase().as_str(),
                    "genome_bed_track" | "genome_bigwig_track" | "blast_hit_track"
                )
            })
    }

    pub fn is_regulatory_feature(feature: &Feature) -> bool {
        let kind = feature.kind.to_string().to_ascii_uppercase();
        Self::is_tfbs_feature(feature)
            || kind.contains("REGULATORY")
            || Self::is_track_feature(feature)
    }

    fn feature_qualifier_text(feature: &Feature, key: &str) -> Option<String> {
        feature
            .qualifier_values(key.into())
            .map(str::trim)
            .find(|value| !value.is_empty())
            .map(ToOwned::to_owned)
    }

    fn first_nonempty_qualifier(feature: &Feature, keys: &[&str]) -> Option<String> {
        for key in keys {
            if let Some(value) = Self::feature_qualifier_text(feature, key) {
                return Some(value);
            }
        }
        None
    }

    pub fn tfbs_group_label(feature: &Feature) -> Option<String> {
        if !Self::is_tfbs_feature(feature) {
            return None;
        }
        let tf_id = Self::feature_qualifier_text(feature, "tf_id");
        let tf_name = Self::first_nonempty_qualifier(
            feature,
            &["bound_moiety", "standard_name", "gene", "name"],
        );
        match (tf_name, tf_id) {
            (Some(name), Some(id)) => {
                if name.eq_ignore_ascii_case(&id) {
                    Some(name)
                } else {
                    Some(format!("{name} ({id})"))
                }
            }
            (Some(name), None) => Some(name),
            (None, Some(id)) => Some(id),
            (None, None) => {
                let fallback = Self::feature_name(feature);
                if fallback.trim().is_empty() {
                    Some("TFBS".to_string())
                } else {
                    Some(fallback)
                }
            }
        }
    }

    pub fn is_restriction_site_feature(feature: &Feature) -> bool {
        let kind = feature.kind.to_string().to_ascii_uppercase();
        if kind == "RESTRICTION_SITE" || kind.contains("RESTRICTION") {
            return true;
        }
        for key in [
            "enzyme",
            "restriction_enzyme",
            "enzyme_name",
            "rebase_enzyme",
        ] {
            if feature.qualifier_values(key.into()).next().is_some() {
                return true;
            }
        }
        false
    }

    pub fn restriction_site_group_label(feature: &Feature) -> Option<String> {
        if !Self::is_restriction_site_feature(feature) {
            return None;
        }
        Self::first_nonempty_qualifier(
            feature,
            &[
                "enzyme",
                "restriction_enzyme",
                "enzyme_name",
                "rebase_enzyme",
                "label",
                "name",
            ],
        )
        .or_else(|| {
            let fallback = Self::feature_name(feature);
            if fallback.trim().is_empty() {
                Some("Restriction site".to_string())
            } else {
                Some(fallback)
            }
        })
    }

    pub fn track_group_label(feature: &Feature) -> Option<String> {
        if !Self::is_track_feature(feature) {
            return None;
        }
        Self::first_nonempty_qualifier(feature, &["gentle_track_name", "name", "label"])
            .or_else(|| Some("Unnamed track".to_string()))
    }

    fn feature_qualifier_f64(feature: &Feature, key: &str) -> Option<f64> {
        feature
            .qualifier_values(key.into())
            .next()
            .and_then(|v| v.trim().parse::<f64>().ok())
    }

    pub fn tfbs_feature_passes_display_filter(
        feature: &Feature,
        criteria: TfbsDisplayCriteria,
    ) -> bool {
        if !Self::is_tfbs_feature(feature) {
            return true;
        }
        if criteria.use_llr_bits {
            let Some(value) = Self::feature_qualifier_f64(feature, "llr_bits") else {
                return false;
            };
            if value < criteria.min_llr_bits {
                return false;
            }
        }
        if criteria.use_llr_quantile {
            let Some(value) = Self::feature_qualifier_f64(feature, "llr_quantile") else {
                return false;
            };
            if value < criteria.min_llr_quantile {
                return false;
            }
        }
        if criteria.use_true_log_odds_bits {
            let value = Self::feature_qualifier_f64(feature, "true_log_odds_bits")
                .or_else(|| Self::feature_qualifier_f64(feature, "log_odds_ratio_bits"));
            let Some(value) = value else {
                return false;
            };
            if value < criteria.min_true_log_odds_bits {
                return false;
            }
        }
        if criteria.use_true_log_odds_quantile {
            let value = Self::feature_qualifier_f64(feature, "true_log_odds_quantile")
                .or_else(|| Self::feature_qualifier_f64(feature, "log_odds_ratio_quantile"));
            let Some(value) = value else {
                return false;
            };
            if value < criteria.min_true_log_odds_quantile {
                return false;
            }
        }
        true
    }

    pub fn is_source_feature(feature: &Feature) -> bool {
        feature.kind.to_string().to_ascii_uppercase() == "SOURCE"
    }

    pub fn feature_name(feature: &Feature) -> String {
        for k in [
            "label",
            "name",
            "standard_name",
            "gene",
            "protein_id",
            "product",
            "region_name",
            "bound_moiety",
        ] {
            let label_text = match feature.qualifier_values(k.into()).next() {
                Some(s) => s.to_owned(),
                None => continue,
            };
            if !label_text.trim().is_empty() {
                return label_text;
            }
        }

        // Keep dense regulatory tracks readable: don't use coordinate-only labels there.
        if Self::is_regulatory_feature(feature) {
            return String::new();
        }

        match feature.location.find_bounds() {
            Ok((from, to)) => format!("{from}..{to}"),
            Err(_) => String::new(),
        }
    }

    pub fn feature_range_text(feature: &Feature) -> String {
        let Ok((from, to_exclusive)) = feature.location.find_bounds() else {
            return String::new();
        };
        if from < 0 || to_exclusive < 0 {
            return String::new();
        }
        let len = to_exclusive.saturating_sub(from);
        let start_1based = from.saturating_add(1);
        format!("{start_1based}..{to_exclusive} ({len} bp)")
    }
}

impl Widget for RenderDna {
    fn ui(self, ui: &mut Ui) -> Response {
        let available = ui.available_size_before_wrap();
        let mut width = if available.x.is_finite() {
            available.x
        } else {
            1.0
        };
        let mut height = if available.y.is_finite() {
            available.y
        } else {
            1.0
        };
        if width <= 0.0 {
            width = ui.max_rect().width().max(1.0);
        }
        if height <= 0.0 {
            height = ui.max_rect().height().max(1.0);
        }
        let safe_size = Vec2::new(width.clamp(1.0, 100_000.0), height.clamp(1.0, 100_000.0));
        let (rect, response) = ui.allocate_exact_size(safe_size, Sense::click());
        let builder = egui::UiBuilder::new()
            .max_rect(rect)
            .layout(ui.layout().clone());
        let mut child_ui = ui.new_child(builder);
        self.render(&mut child_ui, rect);
        response
    }
}
